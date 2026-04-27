#pragma once

#include <amdis/LocalOperators.hpp>

#include <dune/geometry/quadraturerules.hh>
#include <dune/common/referencehelper.hh>

namespace AMDiS {

/**
 * \addtogroup operators
 * @{
 **/

namespace tag {

  template <class Surface, class... GridFunctions>
  struct streamfunction_rhs {

    template <class GF>
    using LocalFunction = std::decay_t<decltype(localFunction(Dune::resolveRef(std::declval<GF const&>())))>;


    int kh = -1;
    Surface const& surface;
    std::tuple<GridFunctions...> gridFcts_;
    std::tuple<LocalFunction<GridFunctions>...> localFcts_;

    streamfunction_rhs(int kh, Surface const& surface, GridFunctions const&... gridFcts) : kh(kh), surface(surface)
      , gridFcts_(gridFcts...)
      , localFcts_(localFunction(Dune::resolveRef(gridFcts))...)
       {}
  };

} // end namespace tag


template <class Surface, class... GridFunctions>
class StreamFunctionRhsOperator
{
  using Self = StreamFunctionRhsOperator;

  template <class GF>
  using LocalFunction = std::decay_t<decltype(localFunction(Dune::resolveRef(std::declval<GF const&>())))>;

public:
  StreamFunctionRhsOperator(tag::streamfunction_rhs<Surface,GridFunctions...> tag) 
    : kh_(tag.kh)
    , surface_(tag.surface)
    , gridFcts_(tag.gridFcts_)
    , localFcts_(tag.localFcts_)
  {}

  // get a local functions
  auto& H()           { return std::get<0>(localFcts_); }
  auto& H() const     { return std::get<0>(localFcts_); }

  // evaluate a local function in a local coordinate x
  auto H(Dune::FieldVector<double,2> const& x) const     { return H()(x); }


  template <class CG, class Node, class Quad, class LocalFct, class Vec>
  void assemble(CG const& contextGeo, Node const& node, Quad const& /*quad*/,
                LocalFct const& localFct, Vec& elementVector) const
  {
    using namespace Dune::Indices;

    static_assert(CG::dim == 2);
    static_assert(CG::dow == 3);

    auto const& phiNode0 = node.child(_0);

    auto const& psiNode0 = node.child(_1);

    auto const& vnNode0 = node.child(_3);

    auto const& HNode0 = node.child(_5);

    std::size_t numPhiLocalFE = phiNode0.finiteElement().size();
    std::size_t numPsiLocalFE = psiNode0.finiteElement().size();
    [[maybe_unused]] std::size_t numVnLocalFE = vnNode0.finiteElement().size();
    std::size_t numHLocalFE = HNode0.finiteElement().size();

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

    // bind the local functions to the current element
    std::apply([&element](auto&... lf) { (lf.bind(element),...); }, localFcts_);
    auto gradH = derivative(H());
    gradH.bind(element);

    using ShapeGradients = typename LFE::Traits::LocalBasisType::Traits::JacobianType;
    using ShapeValues = typename LFE::Traits::LocalBasisType::Traits::RangeType;
    std::vector<ShapeGradients> nShapeGradients1;
    std::vector<ShapeValues> nShapeValues1;
    std::vector<GlobalCoordinate> nGradients1;

    for (std::size_t iq = 0; iq < quad.size(); ++iq)
    {
      auto const& qp = quad[iq];

      // The transposed inverse Jacobian of the map from the reference element to the element
      auto const& jacobian = contextGeo.geometry().jacobianInverseTransposed(qp.position());

      // The multiplicative factor in the integral transformation formula
      auto const dS = contextGeo.geometry().integrationElement(qp.position()) * qp.weight();
      
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
        jacobian.mv(nShapeGradients1[i][0], nGradients1[i]);

      // mean curvature gradient evaluated at QP
      [[maybe_unused]] auto const HAtQP = H(qp.position());
      auto const gradHAtQP = gradH(qp.position());

      // normal gradient evaluated at QP
      FieldMatrix<T,3,3> S(0);
      for (size_t i = 0; i < nGradients1.size(); ++i)
        for (int r = 0; r < 3; ++r)
          for (int s = 0; s < 3; ++s)
            S[r][s] += nGradients1[i][s] * nCoefficients1[i][r];
      S /= nh1_nrm;
      S.rightmultiply(Ph1);

      // Compute the shape function gradients on the real element
      [[maybe_unused]] auto const& phiShapeValues = phiNode0.localBasisValuesAt(qp.position());
      [[maybe_unused]] auto const& psiShapeValues = psiNode0.localBasisValuesAt(qp.position());
      [[maybe_unused]] auto const& vnShapeValues = vnNode0.localBasisValuesAt(qp.position());
      auto const& HShapeValues = HNode0.localBasisValuesAt(qp.position());

      auto const& phiShapeGradients = phiNode0.localBasisJacobiansAt(qp.position());
      auto const& psiShapeGradients = psiNode0.localBasisJacobiansAt(qp.position());
      auto const& vnShapeGradients = vnNode0.localBasisJacobiansAt(qp.position());
      auto const& HShapeGradients = HNode0.localBasisJacobiansAt(qp.position());
      std::vector<GlobalCoordinate> phiGradients(phiShapeGradients.size());
      std::vector<GlobalCoordinate> psiGradients(psiShapeGradients.size());
      std::vector<GlobalCoordinate> vnGradients(vnShapeGradients.size());
      std::vector<GlobalCoordinate> HGradients(HShapeGradients.size());
      for (std::size_t i = 0; i < phiGradients.size(); ++i)
        jacobian.mv(phiShapeGradients[i][0], phiGradients[i]);
      for (std::size_t i = 0; i < psiGradients.size(); ++i)
        jacobian.mv(psiShapeGradients[i][0], psiGradients[i]);
      for (std::size_t i = 0; i < vnGradients.size(); ++i)
        jacobian.mv(vnShapeGradients[i][0],    vnGradients[i]);
      for (std::size_t i = 0; i < HGradients.size(); ++i)
        jacobian.mv(HShapeGradients[i][0],    HGradients[i]);

      auto const exprValue = localFct(qp.position());

      auto H = S[0][0] + S[1][1] + S[2][2];

      // S^2
      FieldMatrix<T,3,3> SS(S);
      SS.leftmultiply(S);
      // tr(S^2)
      auto const trSS = SS[0][0] + SS[1][1] + SS[2][2];
      // K = tr(S)^2 - tr(S^2)
      auto const K = 0.5*(H*H - trSS);

      // Compute the actual vector entries

      // <f_T, curl_S varphi>
      for (std::size_t i = 0; i < numPhiLocalFE; ++i) {
        std::size_t local_i = phiNode0.localIndex(i);
        elementVector[local_i] += dot(exprValue, cross(nh1, phiGradients[i])) * dS;
      }
      // <f_T, grad zeta>
      for (std::size_t i = 0; i < numPsiLocalFE; ++i) {
        std::size_t local_i = psiNode0.localIndex(i);
        elementVector[local_i] += dot(exprValue, psiGradients[i]) * dS;
      }
      // // <f_N, y_N>
      T kappa = Parameters::get<T>("parameters->bending modulus").value_or(1.0);
      // f_N = kappa (Delta H + 0.5*H^3 - 2H K)
      // <f_N, y_N> = -kappa <grad H, grad y_N> + kappa <0.5*H^3 - 2H K, y_N>
      auto const zero_order_term = 0.5*H*H*H - 2*H*K;
      for (std::size_t i = 0; i < numVnLocalFE; ++i) {
        std::size_t local_i = vnNode0.localIndex(i);
        elementVector[local_i] -= kappa * dot(gradHAtQP, vnGradients[i]) * dS;
        elementVector[local_i] += kappa * zero_order_term * vnShapeValues[i] * dS;
      }

      // ---== mean curvature ==---
      // -<grad X, S h> = -<P, Sh>
      T Ph1_S(0);
      for (int r = 0; r < 3; ++r)
        for (int s = 0; s < 3; ++s)
          Ph1_S += Ph1[r][s] * S[r][s];
          
      for (std::size_t i = 0; i < numHLocalFE; ++i) {
        std::size_t local_i = HNode0.localIndex(i);
        elementVector[local_i] -= Ph1_S * HShapeValues[i] * dS;
      }

    }
  }

  private:
  int kh_;
  Surface const& surface_;
  std::tuple<GridFunctions...> gridFcts_;
  mutable std::tuple<LocalFunction<GridFunctions>...> localFcts_;
};

template <class Surface, class... GridFunctions, class LC>
struct GridFunctionOperatorRegistry<tag::streamfunction_rhs<Surface,GridFunctions...>, LC>
{
  static constexpr int degree = 1;
  using type = StreamFunctionRhsOperator<Surface,GridFunctions...>;
};

/**
 * @}
 **/

} // end namespace AMDiS
