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
  struct bgnrhs
  {

    int kh = -1;      // Order of the weingarten-map term
    double dt = 0.0;  // time step
    Surface const& surface;

    bgnrhs(int kh, double dt, Surface const& surface) : kh(kh), dt(dt), surface(surface) {}
  };

} // end namespace tag

template <class Surface>
class BgnRhsOperator
{
  using Self = BgnRhsOperator;
  using ctype = double;

public:
  BgnRhsOperator(tag::bgnrhs<Surface> tag)
    : kh_(tag.kh)
    , dt_(tag.dt)
    , surface_(tag.surface)
  {}

  template <class CG, class Node, class Quad, class LocalFct, class Vec>
  void assemble(CG const& contextGeo, Node const& node,
                Quad const& /*quad*/, LocalFct const& localFct, Vec& elementVector) const
  {
      static_assert(CG::dim == 2);
      static_assert(CG::dow == 3);

      auto const& VecNode = node.child(0_c, 0);
      auto const& localVecFE = VecNode.finiteElement();
      std::size_t sizeVec = localVecFE.size();
      auto const& ScaNode = node.child(1_c);
      auto const& localScaFE = ScaNode.finiteElement();
      std::size_t sizeSca = localScaFE.size();

      // coefficients of normal vector field in lagrange basis
      using T = typename CG::Geometry::ctype;
      using LFE = Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, 2, T, T>;
      using GlobalCoordinate = typename CG::Geometry::GlobalCoordinate;
      std::vector<GlobalCoordinate> nCoefficients1;

      using ShapeValues = typename LFE::Traits::LocalBasisType::Traits::RangeType;
      using ShapeGradients = typename LFE::Traits::LocalBasisType::Traits::JacobianType;
      std::vector<ShapeValues> nShapeValues1;
      std::vector<ShapeGradients> nShapeGradients;
      std::vector<GlobalCoordinate> gradients;

      auto const& element = contextGeo.element();
      auto const& geometry = contextGeo.geometry().impl();
      auto const& quad = Dune::QuadratureRules<T, 2>::rule(element.type(), 20);

      for (std::size_t iq = 0; iq < quad.size(); ++iq)
      {
        auto const& qp = quad[iq];

        // The transposed inverse Jacobian of the map from the reference element to the element
        auto const& jacobian = geometry.jacobianInverseTransposed(qp.position());

        // The multiplicative factor in the integral transformation formula
        auto const dSh = geometry.integrationElement(qp.position()) * qp.weight();

        // normal vector evaluated at QP
        GlobalCoordinate nh = geometry.normal(qp.position());

        FieldMatrix<T,3,3> Ph;
        for (int r = 0; r < 3; ++r)
          for (int s = 0; s < 3; ++s)
            Ph[r][s] = (r == s ? 1 : 0) - nh[r]*nh[s];

        auto vnAtQP = 1.0;
        auto const exprValue = localFct(qp.position());

        // ----- evaluation of values and gradients of trial/test basis functions ------------

        [[maybe_unused]] auto const& vecShapeValues = VecNode.localBasisValuesAt(qp.position());
        auto const& scaShapeValues = ScaNode.localBasisValuesAt(qp.position());
        auto const& shapeGradients = VecNode.localBasisJacobiansAt(qp.position());
        std::vector<GlobalCoordinate> gradients(shapeGradients.size());
        for (std::size_t i = 0; i < gradients.size(); ++i)
          jacobian.mv(shapeGradients[i][0], gradients[i]);  

        // ---== <Y * n, h> = <vn*dt, h> ==---
        // <vn*dt, h>
        for (std::size_t i = 0; i < sizeSca; ++i) {
            auto const local_i = node.child(1_c).localIndex(i);
            elementVector[local_i] += vnAtQP * dt_ * scaShapeValues[i] * dSh;
        }

        // ---== <H*n, y> + <grad Y, grad y> = -<grad X, grad y> ==---

        // -<grad X, grad y> = -<P, grad y>
        for (std::size_t i = 0; i < sizeVec; ++i) {
          for (int r = 0; r < 3; ++r) {
            for (int s = 0; s < 3; ++s) {
              std::size_t row = node.child(0_c,s).localIndex(i);

              // elementVector[row] -= exprValue[s][r] * gradients[i][r] * dSh;
              elementVector[row] -= Ph[s][r] * gradients[i][r] * dSh;
            }
          }
        }
      }
  }

private:
  int kh_;
  double dt_;
  Surface const& surface_;
};


template <class Surface, class LC>
struct GridFunctionOperatorRegistry<tag::bgnrhs<Surface>, LC>
{
  static constexpr int degree = 2;
  using type = BgnRhsOperator<Surface>;
};


/**
 * @}
 **/

} // end namespace AMDiS
