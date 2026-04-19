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

template <class Surface, class VnGridFct, class XGridFct>
class BgnRhsOperator
{
  using Self = BgnRhsOperator;

  template <class VnLocalFct, class XLocalFct>
  class BgnRhsLocalOperator
  {
    using ctype = double;
    using LFECache = Dune::LagrangeLFECache<ctype, ctype, 2>;
    using LFE = typename LFECache::FiniteElementType;

  public:
    BgnRhsLocalOperator(int kh, double dt, Surface const& surface
                        , VnLocalFct const& vnLocalFct
                        , XLocalFct const& xLocalFct)
      : lfeCache_kh_(std::max(1,kh-1))
      , dt_(dt)
      , surface_(surface)
      , vnLocalFct_(vnLocalFct)
      , xLocalFct_(xLocalFct)
    {}

    template <class Element>
    void bind(Element const& element){
      vnLocalFct_.bind(element);
      xLocalFct_.bind(element);
    }

    void unbind(){
      vnLocalFct_.unbind();
      xLocalFct_.unbind();
    }

    template <class CG, class Node, class Vec>
    void assemble(CG const& contextGeo, Node const& node, Vec& elementVector) const
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

        auto vnAtQP = vnLocalFct_(qp.position());

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

        // -<grad X, grad y> =? -<I, grad y>
        for (std::size_t i = 0; i < sizeVec; ++i) {
          for (int r = 0; r < 3; ++r) {
            std::size_t row = node.child(0_c,r).localIndex(i);

            elementVector[row] -= gradients[i][r] * dSh;
          }
        }
      }
    }

  private:
    mutable LFECache lfeCache_kh_;
    double dt_;
    Surface const& surface_;
    VnLocalFct vnLocalFct_;
    XLocalFct xLocalFct_;
  };

public:
  BgnRhsOperator(int kh, double dt, Surface const& surface, VnGridFct const& vnGridFct, XGridFct const& xGridFct)
    : kh_(kh)
    , dt_(dt)
    , surface_(surface)
    , vnGridFct_(vnGridFct)
    , xGridFct_(xGridFct)
  {}

  template <class GridView>
  void update(GridView const&) { /* do nothing */ }

  friend auto localOperator(BgnRhsOperator const& op)
  {
    return BgnRhsLocalOperator{op.kh_, op.dt_, op.surface_,
      localFunction(op.vnGridFct_),
      localFunction(op.xGridFct_)};
  }

  int kh_;
  double dt_;
  Surface const& surface_;
  VnGridFct vnGridFct_;
  xGridFct xGridFct_;
};


template <class Surface, class Vn, class X>
struct PreBgnRhsOperator
{
  int kh;
  double dt;
  Surface const& surface;
  Vn vn;
  X x;
};

template <class Surface, class Vn, class X>
auto bgnrhs(int kh, double dt, Surface const& surface, Vn const& vn, X const& x)
{
  return PreBgnRhsOperator<Surface, Vn, X>{kh, dt, surface, vn, x};
}

template <class Context, class... T, class GridView>
auto makeOperator(PreBgnRhsOperator<T...> pre, GridView const& gridView)
{
  return BgnRhsOperator{pre.kh, pre.dt, pre.surface, 
    makeGridFunction(std::move(pre.vn), gridView),
    makeGridFunction(std::move(pre.x), gridView)};
}

/**
 * @}
 **/

} // end namespace AMDiS
