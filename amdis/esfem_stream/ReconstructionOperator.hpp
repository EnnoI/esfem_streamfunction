#pragma once

#include <amdis/LocalOperators.hpp>
#include <amdis/common/Literals.hpp>

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/lagrange/lfecache.hh>

namespace AMDiS {

/**
 * \addtogroup operators
 * @{
 **/

template <class Surface, class Data, class PhiGridFct, class PsiGridFct, class FGridFct>
class ReconstructionOperator
{
  using Self = ReconstructionOperator;

  template <class PhiLocalFct, class PsiLocalFct, class FLocalFct>
  class ReconstructionLocalOperator
  {
    using ctype = double;
    using LFECache = Dune::LagrangeLFECache<ctype, ctype, 2>;
    using LFE = typename LFECache::FiniteElementType;
    using GradPhiLocalFct = TYPEOF(derivativeOf(std::declval<PhiLocalFct>(), tag::gradient{}));
    using GradPsiLocalFct = TYPEOF(derivativeOf(std::declval<PsiLocalFct>(), tag::gradient{}));

  public:
    ReconstructionLocalOperator(int kh, Surface const& surface, Data& data,
                                PhiLocalFct const& phiLocalFct, PsiLocalFct const& psiLocalFct, FLocalFct const& fLocalFct)
      : lfeCache_kh_(std::max(1,kh-1))
      , surface_(surface)
      , data_(data)
      , gradPhiLocalFct_(derivativeOf(phiLocalFct, tag::gradient{}))
      , gradPsiLocalFct_(derivativeOf(psiLocalFct, tag::gradient{}))
      , fLocalFct_(fLocalFct)
    {}

    template <class Element>
    void bind(Element const& element)
    {
      gradPhiLocalFct_.bind(element);
      gradPsiLocalFct_.bind(element);
      fLocalFct_.bind(element);
    }

    void unbind()
    {
      gradPhiLocalFct_.unbind();
      gradPsiLocalFct_.unbind();
      fLocalFct_.unbind();
    }

    template <class CG, class Node, class Vec>
    void assemble(CG const& contextGeo, Node const& node, Vec& elementVector) const
    {
      static_assert(CG::dim == 2);
      static_assert(CG::dow == 3);

      auto const& localVelFE = node.child(0_c,0).finiteElement();
      std::size_t sizeVel = localVelFE.size();

      // coefficients of normal vector field in lagrange basis
      using T = typename CG::Geometry::ctype;
      using GlobalCoordinate = typename CG::Geometry::GlobalCoordinate;
      using ShapeGradients = typename LFE::Traits::LocalBasisType::Traits::JacobianType;
      using ShapeValues = typename LFE::Traits::LocalBasisType::Traits::RangeType;

      auto const& element = contextGeo.element();
      auto const& geometry = contextGeo.geometry().impl();
      auto const& quad = Dune::QuadratureRules<T, 2>::rule(element.type(), 20);

      auto const& localFE_kh = lfeCache_kh_.get(contextGeo.type());
      std::vector<GlobalCoordinate> nCoefficients1;

      // Interpolate discrete normal into FE space
      localFE_kh.localInterpolation().interpolate(
        [&](auto const& local) -> GlobalCoordinate
        {
          return geometry.normal(local);
        },
        nCoefficients1
        );

      std::vector<ShapeGradients> nShapeGradients1;
      std::vector<ShapeValues> nShapeValues1;
      std::vector<GlobalCoordinate> nGradients1;

      for (std::size_t iq = 0; iq < quad.size(); ++iq)
      {
        auto const& qp = quad[iq];

        // The transposed inverse Jacobian of the map from the reference element to the element
        [[maybe_unused]] auto const& jacobian = geometry.jacobianInverseTransposed(qp.position());

        // The multiplicative factor in the integral transformation formula
        auto const dSh = geometry.integrationElement(qp.position()) * qp.weight();

        // ------- evaluation of normal vector and its derivatives at quadrature point  ------
        localFE_kh.localBasis().evaluateFunction(qp.position(), nShapeValues1);
        localFE_kh.localBasis().evaluateJacobian(qp.position(), nShapeGradients1);

        // normal vector evaluated at QP
        GlobalCoordinate nh = geometry.normal(qp.position());
        auto nh_nrm = nh.two_norm();

        FieldMatrix<T,3,3> Ph;
        for (int r = 0; r < 3; ++r)
          for (int s = 0; s < 3; ++s)
            Ph[r][s] = (r == s ? 1 : 0) - nh[r]*nh[s];

        // Compute the shape function gradients on the real element
        nGradients1.resize(nShapeGradients1.size());
        for (size_t i = 0; i < nGradients1.size(); ++i)
          jacobian.mv(nShapeGradients1[i][0], nGradients1[i]);

        // normal gradient evaluated at QP
        [[maybe_unused]] FieldMatrix<T,3,3> S(0);
        for (size_t i = 0; i < nGradients1.size(); ++i)
          for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
              S[r][s] += nGradients1[i][s] * nCoefficients1[i][r];
        S /= nh_nrm;
        S.leftmultiply(Ph);
        S.rightmultiply(Ph);

        // H = tr(S)
        // auto const H = S[0][0] + S[1][1] + S[2][2];

        [[maybe_unused]] auto fAtQP = fLocalFct_(qp.position());
        auto gradPhiAtQP = gradPhiLocalFct_(qp.position());
        auto gradPsiAtQP = Ph * gradPsiLocalFct_(qp.position());
        auto curlPhiAtQP = Ph * cross(nh, gradPhiAtQP);

        // ----- evaluation of values and gradients of trial/test basis functions ------------

        auto const& shapeValues = node.child(0_c,0).localBasisValuesAt(qp.position());

        for (std::size_t j = 0; j < sizeVel; ++j) {
          for (int i = 0; i < 3; ++i) {
            std::size_t row = node.child(0_c,i).localIndex(j);

            // (curl(phi) + grad(psi), vh)
            elementVector[row] += (curlPhiAtQP[i] + gradPsiAtQP[i]) * shapeValues[j] * dSh;
          }
        }

        // TODO: pressure reconstruction
      }
    }

  private:
    mutable LFECache lfeCache_kh_;

    Surface const& surface_;
    Data& data_;

    GradPhiLocalFct gradPhiLocalFct_;
    GradPsiLocalFct gradPsiLocalFct_;
    FLocalFct fLocalFct_;
  };

public:
  ReconstructionOperator(int kh, Surface const& surface, Data& data,
                         PhiGridFct const& phiGridFct, PsiGridFct const& psiGridFct, FGridFct const& fGridFct)
    : kh_(kh)
    , surface_(surface)
    , data_(data)
    , phiGridFct_(phiGridFct)
    , psiGridFct_(psiGridFct)
    , fGridFct_(fGridFct)
  {}

  template <class GridView>
  void update(GridView const&) { /* do nothing */ }

  friend auto localOperator(ReconstructionOperator const& op)
  {
    return ReconstructionLocalOperator{op.kh_, op.surface_, op.data_,
      localFunction(op.phiGridFct_),
      localFunction(op.psiGridFct_),
      localFunction(op.fGridFct_)};
  }

  int kh_;

  Surface const& surface_;
  Data& data_;

  PhiGridFct phiGridFct_;
  PsiGridFct psiGridFct_;
  FGridFct fGridFct_;
};


template <class Surface, class Data, class Phi, class Psi, class F>
struct PreReconstructionOperator
{
  int kh;
  Surface const& surface;
  Data& data;
  Phi phi;
  Psi psi;
  F f;
};

template <class Surface, class Data, class Phi, class Psi, class F>
auto reconstruction(int kh, Surface const& surface, Data& data, Phi const& phi, Psi const& psi, F const& f)
{
  return PreReconstructionOperator<Surface,Data,Phi,Psi,F>{kh, surface, data, phi, psi, f};
}

template <class Context, class... T, class GridView>
auto makeOperator(PreReconstructionOperator<T...> pre, GridView const& gridView)
{
  return ReconstructionOperator{pre.kh, pre.surface, pre.data,
    makeGridFunction(std::move(pre.phi), gridView),
    makeGridFunction(std::move(pre.psi), gridView),
    makeGridFunction(std::move(pre.f), gridView)};
}

/**
 * @}
 **/

} // end namespace AMDiS
