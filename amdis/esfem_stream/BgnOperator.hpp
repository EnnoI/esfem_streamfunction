#pragma once

#include <amdis/LocalOperators.hpp>

#include <dune/geometry/quadraturerules.hh>
// #include <dune/common/referencehelper.hh>

namespace AMDiS {

/**
 * \addtogroup operators
 * @{
 **/

namespace tag {

  template <class Surface>
  struct bgn
  {

    int kh = -1;      // Order of the weingarten-map term
    Surface const& surface;

    bgn(int kh, Surface const& surface) : kh(kh), surface(surface) {}
  };

} // end namespace tag

template <class Surface>
class BgnOperator
{
  using Self = BgnOperator;
  using ctype = double;

public:
  BgnOperator(tag::bgn<Surface> tag)
    : kh_(tag.kh)
    , surface_(tag.surface)
  {}

  template <class CG, class Node, class Quad, class LocalFct, class Mat>
  void assemble(CG const& contextGeo, Node const& node, Node const& colNode,
                Quad const& /*quad*/, LocalFct const& /*localFct*/, Mat& elementMatrix) const
  {
    using namespace Dune::Indices;
    static_assert(CG::dim == 2);
    static_assert(CG::dow == 3);

    auto const& yNode0 = node.child(_0);
    auto const& yNode1 = colNode.child(_0);

    auto const& hNode0 = node.child(_1);
    auto const& hNode1 = colNode.child(_1);

    std::size_t numYFE = yNode0.child(0).finiteElement().size();
    std::size_t numHFE = hNode0.finiteElement().size();

    using GlobalCoordinate = typename CG::Geometry::GlobalCoordinate;
    using T = typename CG::Geometry::ctype;

    auto const& element = contextGeo.element();
    auto const& geometry = contextGeo.geometry().impl();
    auto const& quad = Dune::QuadratureRules<T, 2>::rule(element.type(), 20);

    using LFE = Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, 2, T, T>;

    auto localFE_kh = LFE{element.type(), (unsigned int)(kh_)};
    std::vector<GlobalCoordinate> nCoefficients1;

    // Interpolate discrete normal into FE space
    localFE_kh.localInterpolation().interpolate(
      [&](auto const& local) -> GlobalCoordinate
      {
        return geometry.normal(local);
      },
      nCoefficients1
    );

    using ShapeGradients = typename LFE::Traits::LocalBasisType::Traits::JacobianType;
    using ShapeValues = typename LFE::Traits::LocalBasisType::Traits::RangeType;
    std::vector<ShapeGradients> nShapeGradients1;
    std::vector<ShapeValues> nShapeValues1;
    std::vector<GlobalCoordinate> nGradients1;

    for (std::size_t iq = 0; iq < quad.size(); ++iq)
    {
      auto const& qp = quad[iq];
      auto const dSh = geometry.integrationElement(qp.position()) * qp.weight();
      auto const& Jh = geometry.jacobianInverseTransposed(qp.position());

      // ------- evaluation of normal vector and its derivatives at quadrature point  ------
      localFE_kh.localBasis().evaluateFunction(qp.position(), nShapeValues1);
      localFE_kh.localBasis().evaluateJacobian(qp.position(), nShapeGradients1);

      GlobalCoordinate nh1(0);
      for (std::size_t i = 0; i < nShapeValues1.size(); ++i)
        nh1.axpy(nShapeValues1[i], nCoefficients1[i]);

      auto nh1_nrm = nh1.two_norm();
      nh1 /= nh1_nrm;

      FieldMatrix<T,3,3> Ph1;
      for (int r = 0; r < 3; ++r) {
        for (int s = 0; s < 3; ++s) {
          Ph1[r][s] = (r == s ? 1 : 0) - nh1[r]*nh1[s];
        }
      }

      // Compute the shape function gradients on the real element
      nGradients1.resize(nShapeGradients1.size());
      for (size_t i = 0; i < nGradients1.size(); ++i)
        Jh.mv(nShapeGradients1[i][0], nGradients1[i]);

      // normal gradient evaluated at QP
      FieldMatrix<T,3,3> S(0);
      for (size_t i = 0; i < nGradients1.size(); ++i)
        for (int r = 0; r < 3; ++r)
          for (int s = 0; s < 3; ++s)
            S[r][s] += nGradients1[i][s] * nCoefficients1[i][r];
      S /= nh1_nrm;
      S.rightmultiply(Ph1);

      // H = tr(S)
      [[maybe_unused]] auto H = S[0][0] + S[1][1] + S[2][2];

      // ----- evaluation of values and gradients of trial/test basis functions ------------

      auto const& yShapeValues = yNode0.child(0).localBasisValuesAt(qp.position());
      auto const& hShapeValues = hNode0.localBasisValuesAt(qp.position());
      auto const& yShapeGradients = yNode0.child(0).localBasisJacobiansAt(qp.position());
      auto const& hShapeGradients = hNode0.localBasisJacobiansAt(qp.position()); 
      std::vector<GlobalCoordinate> yGradients(yShapeGradients.size());
      std::vector<GlobalCoordinate> hGradients(hShapeGradients.size());
      for (std::size_t i = 0; i < yGradients.size(); ++i)
        Jh.mv(yShapeGradients[i][0], yGradients[i]);  
      for (std::size_t i = 0; i < hGradients.size(); ++i)
        Jh.mv(hShapeGradients[i][0], hGradients[i]);  

      // ----- compute the actual matrix entries --------------

      // ---== <Y * n, h> = <vn*dt, h> ==---
      // <Y * n, h>
      for (std::size_t i = 0; i < numHFE; ++i) {
        for (std::size_t j = 0; j < numYFE; ++j) {
          for (int r = 0; r < 3; ++r) {
            auto const local_i = hNode0.localIndex(i);
            auto const local_j = yNode1.child(r).localIndex(j);
            elementMatrix[local_i][local_j] += yShapeValues[j] * nh1[r] * hShapeValues[i] * dSh;
          }
        }
      }


      // ---== <H*n, y> + <grad Y, grad y> = -<grad X, grad y> ==---

      // <grad Y, grad y> = tr(grad Y^T grad y)
      for (std::size_t i = 0; i < numYFE; ++i) {
        for (std::size_t j = 0; j < numYFE; ++j) {
          for (int r = 0; r < 3; ++r) {
            for (int s = 0; s < 3; ++s) {
              auto const local_i = yNode0.child(r).localIndex(i);
              auto const local_j = yNode1.child(r).localIndex(j);
              elementMatrix[local_i][local_j] += yGradients[i][s] * yGradients[j][s] * dSh;
            }
          }
        }
      }

      // <H*n, y>
      for (std::size_t i = 0; i < numYFE; ++i) {
        for (std::size_t j = 0; j < numHFE; ++j) {
          for (int r = 0; r < 3; ++r) {
            auto const local_i = yNode0.child(r).localIndex(i);
            auto const local_j = hNode1.localIndex(j);
            elementMatrix[local_i][local_j] += hShapeValues[j] * yShapeValues[i]*nh1[r] * dSh;
          }
        }
      }

    }





  }


private:
  int kh_;
  Surface const& surface_;
};


template <class Surface, class LC>
struct GridFunctionOperatorRegistry<tag::bgn<Surface>, LC>
{
  static constexpr int degree = 2;
  using type = BgnOperator<Surface>;
};


/**
 * @}
 **/

} // end namespace AMDiS
