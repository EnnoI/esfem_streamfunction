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

  template <class Surface, class Data, class... GridFunctions>
  struct stream_function
  {

    template <class GF>
    using LocalFunction = std::decay_t<decltype(localFunction(Dune::resolveRef(std::declval<GF const&>())))>;

    int kg = -1;      // Order of the geometry approximation
    int kh = -1;      // Order of the weingarten-map term
    Surface const& surface;
    Data& data;
    std::tuple<GridFunctions...> gridFcts_;
    std::tuple<LocalFunction<GridFunctions>...> localFcts_;


    stream_function(int kg, int kh, Surface const& surface, Data& data, GridFunctions const&... gridFcts)
      : kg(kg), kh(kh), surface(surface), data(data)
      , gridFcts_(gridFcts...)
      , localFcts_(localFunction(Dune::resolveRef(gridFcts))...)
    {}
  };

} // end namespace tag


template <class Surface, class Data, class... GridFunctions>
class StreamFunctionOperator
{
  using Self = StreamFunctionOperator;
  using ctype = double;
  template <class GF>
  using LocalFunction = std::decay_t<decltype(localFunction(Dune::resolveRef(std::declval<GF const&>())))>;

public:
  StreamFunctionOperator(tag::stream_function<Surface,Data,GridFunctions...> tag)
    : kg_(tag.kg)
    , kh_(tag.kh)
    , surface_(tag.surface)
    , data_(tag.data)
    , gridFcts_(tag.gridFcts_)
    , localFcts_(tag.localFcts_)
  {}

  // get a local functions
  auto& H()           { return std::get<0>(localFcts_); }
  auto& H() const     { return std::get<0>(localFcts_); }

  // evaluate a local function in a local coordinate x
  auto H(Dune::FieldVector<double,2> const& x) const     { return H()(x); }

  template <class CG, class Node, class Quad, class LocalFct, class Mat>
  void assemble(CG const& contextGeo, Node const& node, Node const& colNode,
                Quad const& /*quad*/, LocalFct const& /*localFct*/, Mat& elementMatrix) const
  {
    using namespace Dune::Indices;
    static_assert(CG::dim == 2);
    static_assert(CG::dow == 3);

    auto const& phiNode0 = node.child(_0);
    auto const& phiNode1 = colNode.child(_0);

    auto const& psiNode0 = node.child(_1);
    auto const& psiNode1 = colNode.child(_1);

    auto const& omegaNode0 = node.child(_2);
    auto const& omegaNode1 = colNode.child(_2);

    auto const& vnNode0 = node.child(_3);
    auto const& vnNode1 = colNode.child(_3);

    auto const& pnNode0 = node.child(_4);
    auto const& pnNode1 = colNode.child(_4);

    auto const& HNode0 = node.child(_5);
    auto const& HNode1 = colNode.child(_5);

    std::size_t numPhiLocalFE = phiNode0.finiteElement().size();
    std::size_t numPsiLocalFE = psiNode0.finiteElement().size();
    std::size_t numOmegaLocalFE = omegaNode0.finiteElement().size();
    std::size_t numVnLocalFE = vnNode0.finiteElement().size();
    std::size_t numPnLocalFE = pnNode0.finiteElement().size();
    std::size_t numHLocalFE = HNode0.finiteElement().size();

    using GlobalCoordinate = typename CG::Geometry::GlobalCoordinate;
    using T = ctype;

    auto const& element = contextGeo.element();
    auto const& geometry = contextGeo.geometry().impl();
    auto const& quad = Dune::QuadratureRules<T, 2>::rule(element.type(), 14);

    using LFE = Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, 2, T, T>;

    auto localFE_kh = LFE{element.type(), (unsigned int)(kh_)};
    std::vector<GlobalCoordinate> nCoefficients1;

    // Interpolate discrete normal into FE space
    localFE_kh.localInterpolation().interpolate(
      [&](auto const& local) -> GlobalCoordinate
      {
        // auto const& JT = geometry.jacobianTransposed(local);

        // GlobalCoordinate n = cross(JT[0], JT[1]);

        // n /= n.two_norm();
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

      // H = tr(S)
      auto H = S[0][0] + S[1][1] + S[2][2];
      // S^2
      FieldMatrix<T,3,3> SS(S);
      SS.leftmultiply(S);
      // tr(S^2)
      auto const trSS = SS[0][0] + SS[1][1] + SS[2][2];
      // K = tr(S)^2 - tr(S^2)
      auto const K = 0.5*(H*H - trSS);
      // S : S
      T SdotS(0);
      for (int r = 0; r < 3; ++r)
        for (int s = 0; s < 3; ++s)
          SdotS += S[r][s]*S[r][s];


      // ----- evaluation of values and gradients of trial/test basis functions ------------

      [[maybe_unused]] auto const& phiShapeValues = phiNode0.localBasisValuesAt(qp.position());
      [[maybe_unused]] auto const& psiShapeValues = psiNode0.localBasisValuesAt(qp.position());
      auto const& omegaShapeValues = omegaNode0.localBasisValuesAt(qp.position());
      auto const& vnShapeValues = vnNode0.localBasisValuesAt(qp.position());
      auto const& pnShapeValues = pnNode0.localBasisValuesAt(qp.position());
      auto const& HShapeValues = HNode0.localBasisValuesAt(qp.position());

      auto const& phiShapeGradients = phiNode0.localBasisJacobiansAt(qp.position());
      auto const& psiShapeGradients = psiNode0.localBasisJacobiansAt(qp.position());
      auto const& omegaShapeGradients = omegaNode0.localBasisJacobiansAt(qp.position());
      auto const& vnShapeGradients = vnNode0.localBasisJacobiansAt(qp.position());
      auto const& pnShapeGradients = pnNode0.localBasisJacobiansAt(qp.position());
      auto const& HShapeGradients = HNode0.localBasisJacobiansAt(qp.position());
      std::vector<GlobalCoordinate> phiGradients(phiShapeGradients.size());
      std::vector<GlobalCoordinate> psiGradients(psiShapeGradients.size());
      std::vector<GlobalCoordinate> omegaGradients(omegaShapeGradients.size());
      std::vector<GlobalCoordinate> vnGradients(vnShapeGradients.size());
      std::vector<GlobalCoordinate> pnGradients(pnShapeGradients.size());
      std::vector<GlobalCoordinate> HGradients(HShapeGradients.size());
      for (std::size_t i = 0; i < phiGradients.size(); ++i)
        Jh.mv(phiShapeGradients[i][0], phiGradients[i]);
      for (std::size_t i = 0; i < psiGradients.size(); ++i)
        Jh.mv(psiShapeGradients[i][0], psiGradients[i]);
      for (std::size_t i = 0; i < omegaGradients.size(); ++i)
        Jh.mv(omegaShapeGradients[i][0], omegaGradients[i]);
      for (std::size_t i = 0; i < vnGradients.size(); ++i)
        Jh.mv(vnShapeGradients[i][0],    vnGradients[i]);
      for (std::size_t i = 0; i < pnGradients.size(); ++i)
        Jh.mv(pnShapeGradients[i][0],    pnGradients[i]);  
      for (std::size_t i = 0; i < HGradients.size(); ++i)
        Jh.mv(HShapeGradients[i][0],     HGradients[i]);

      // ----- compute the actual matrix entries --------------

      T mu = Parameters::get<T>("parameters->mu").value_or(1.0);
      T reg = Parameters::get<T>("parameters->kill regularization").value_or(0.0);

      // ---== curl MB ==---
      // -mu<grad w, grad varphi>
      for (std::size_t i = 0; i < numPhiLocalFE; ++i) {
        for (std::size_t j = 0; j < numOmegaLocalFE; ++j) {
          auto const local_i = phiNode0.localIndex(i);
          auto const local_j = omegaNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= mu * dot(omegaGradients[j], phiGradients[i]) * dSh;
        }
      }
      // -<2muK grad phi, grad varphi> + reg * <grad phi, grad varphi>
      for (std::size_t i = 0; i < numPhiLocalFE; ++i) {
        for (std::size_t j = 0; j < numPhiLocalFE; ++j) {
          auto const local_i = phiNode0.localIndex(i);
          auto const local_j = phiNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= (2*mu*K - reg) * dot(phiGradients[j], phiGradients[i]) * dSh; // reg or mu*reg, + or - ?
        }
      }
      // -<2mu P div(v_N S), curl varphi> = -2mu<grad v_N S + v_N grad H, curl varphi>
      // grad y_N^T S = grad v_N^T S
      decltype(vnGradients) grad_y_N_S(vnGradients.size());
      for (std::size_t i = 0; i < grad_y_N_S.size(); ++i) {
        for (std::size_t r = 0; r < 3; ++r)
          for (std::size_t s = 0; s < 3; ++s)
            grad_y_N_S[i][r] += vnGradients[i][s] * S[s][r];
      } 
      decltype(vnGradients) div_y_N_S(vnGradients.size());
      for (std::size_t i = 0; i < div_y_N_S.size(); ++i) {
          div_y_N_S[i] = grad_y_N_S[i] + vnShapeValues[i] * gradHAtQP;
      }
      // -2mu<grad v_N S + v_N grad H, curl varphi>
      for (std::size_t i = 0; i < numPhiLocalFE; ++i) {
        for (std::size_t j = 0; j < numVnLocalFE; ++j) {
          auto const local_i = phiNode0.localIndex(i);
          auto const local_j = vnNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= 2 * mu * dot(Ph1 * div_y_N_S[j], cross(nh1,phiGradients[i])) * dSh;
        }
      }
      // ---== div MB ==---
      // -<2muK grad psi, grad zeta> + reg * <grad psi, grad zeta>
      for (std::size_t i = 0; i < numPsiLocalFE; ++i) {
        for (std::size_t j = 0; j < numPsiLocalFE; ++j) {
          auto const local_i = psiNode0.localIndex(i);
          auto const local_j = psiNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= (2*mu*K - reg) * dot(psiGradients[j], psiGradients[i]) * dSh; // reg or mu*reg, + or - ?
        }
      }
      // <2mu grad(v_N H), grad zeta> = 2mu <grad v_N H + v_N grad H, grad zeta>
      for (std::size_t i = 0; i < numPsiLocalFE; ++i) {
        for (std::size_t j = 0; j < numVnLocalFE; ++j) {
          auto const local_i = psiNode0.localIndex(i);
          auto const local_j = vnNode1.localIndex(j);
          elementMatrix[local_i][local_j] += 2 * mu * dot(vnGradients[j] * HAtQP + vnShapeValues[j] * gradHAtQP, psiGradients[i]) * dSh;
        }
      }
      // -<2 mu P div(v_N S), grad zeta> = -2mu <grad v_N S + v_N grad H, grad zeta>
      for (std::size_t i = 0; i < numPsiLocalFE; ++i) {
        for (std::size_t j = 0; j < numVnLocalFE; ++j) {
          auto const local_i = psiNode0.localIndex(i);
          auto const local_j = vnNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= 2 * mu * dot(Ph1 * div_y_N_S[j], psiGradients[i]) * dSh;
        }
      }
      // <grad p, grad zeta>
      for (std::size_t i = 0; i < numPsiLocalFE; ++i) {
        for (std::size_t j = 0; j < numPnLocalFE; ++j) {
          auto const local_i = psiNode0.localIndex(i);
          auto const local_j = pnNode1.localIndex(j);
          elementMatrix[local_i][local_j] += dot(pnGradients[j], psiGradients[i]) * dSh;
        }
      }

      // ---== normal MB ==---
      // <mu E(u_T), grad(y_N n)> = -2mu<grad psi + curl phi, grad y_N S + y_N grad H>
      // = -2mu<grad psi, grad y_N S + y_N grad H> 
      for (std::size_t i = 0; i < numVnLocalFE; ++i) {
        for (std::size_t j = 0; j < numPsiLocalFE; ++j) {
          auto const local_i = vnNode0.localIndex(i);
          auto const local_j = psiNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= 2 * mu * dot(psiGradients[j], div_y_N_S[i]) * dSh; // + or - ?? - should be correct but + works...
        }
      }
      // - 2mu<curl phi, grad y_N S + y_N grad H>
      for (std::size_t i = 0; i < numVnLocalFE; ++i) {
        for (std::size_t j = 0; j < numPhiLocalFE; ++j) {
          auto const local_i = vnNode0.localIndex(i);
          auto const local_j = phiNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= 2 * mu * dot(cross(nh1,phiGradients[j]), div_y_N_S[i]) * dSh; // + or - ??
        }
      }
      // <2 mu v_N S, y_N S> - reg * <v_N, y_N>
      for (std::size_t i = 0; i < numVnLocalFE; ++i) {
        for (std::size_t j = 0; j < numVnLocalFE; ++j) {
          auto const local_i = vnNode0.localIndex(i);
          auto const local_j = vnNode1.localIndex(j);
          elementMatrix[local_i][local_j] += (2 * mu * SdotS - reg) * vnShapeValues[i] * vnShapeValues[j] * dSh;
        }
      }
      // <-pH, y_N>
      for (std::size_t i = 0; i < numVnLocalFE; ++i) {
        for (std::size_t j = 0; j < numPnLocalFE; ++j) {
          auto const local_i = vnNode0.localIndex(i);
          auto const local_j = pnNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= H * vnShapeValues[i] * pnShapeValues[j] * dSh;
        }
      }

      // ---== vorticity == ---
      // <grad phi, grad delta>
      for (std::size_t i = 0; i < numOmegaLocalFE; ++i) {
        for (std::size_t j = 0; j < numPhiLocalFE; ++j) {
          auto const local_i = omegaNode0.localIndex(i);
          auto const local_j = phiNode1.localIndex(j);
          elementMatrix[local_i][local_j] += dot(phiGradients[j], omegaGradients[i]) * dSh;
        }
      }
      // <omega, delta>
      for (std::size_t i = 0; i < numOmegaLocalFE; ++i) {
        for (std::size_t j = 0; j < numOmegaLocalFE; ++j) {
          auto const local_i = omegaNode0.localIndex(i);
          auto const local_j = omegaNode1.localIndex(j);
          elementMatrix[local_i][local_j] += omegaShapeValues[j] * omegaShapeValues[i] * dSh;
        }
      }

      // ---== incompressibility ==---
      // <grad psi, grad q>
      for (std::size_t i = 0; i < numPnLocalFE; ++i) {
        for (std::size_t j = 0; j < numPsiLocalFE; ++j) {
          auto const local_i = pnNode0.localIndex(i);
          auto const local_j = psiNode1.localIndex(j);
          elementMatrix[local_i][local_j] += dot(psiGradients[j], pnGradients[i]) * dSh;
        }
      }
      // <-v_N H, q>
      for (std::size_t i = 0; i < numPnLocalFE; ++i) {
        for (std::size_t j = 0; j < numVnLocalFE; ++j) {
          auto const local_i = pnNode0.localIndex(i);
          auto const local_j = vnNode1.localIndex(j);
          elementMatrix[local_i][local_j] -= H * vnShapeValues[j] * pnShapeValues[i] * dSh;
        }
      }
      // ---== mean curvature ==---
      // <H, h>
      for (std::size_t i = 0; i < numHLocalFE; ++i) {
        for (std::size_t j = 0; j < numHLocalFE; ++j) {
          auto const local_i = HNode0.localIndex(i);
          auto const local_j = HNode1.localIndex(j);
          elementMatrix[local_i][local_j] += HShapeValues[j] * HShapeValues[i] * dSh;
        }
      }
    }
  }


private:
  int kg_;
  int kh_;
  Surface const& surface_;
  Data& data_;
  std::tuple<GridFunctions...> gridFcts_;
  mutable std::tuple<LocalFunction<GridFunctions>...> localFcts_;
};


template <class Surface, class Data, class... GridFunctions, class LC>
struct GridFunctionOperatorRegistry<tag::stream_function<Surface,Data,GridFunctions...>, LC>
{
  static constexpr int degree = 2;
  using type = StreamFunctionOperator<Surface,Data,GridFunctions...>;
};


/**
 * @}
 **/

} // end namespace AMDiS
