/*
-----
RHS of shape operator computation
-----
*/

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

  template <class... GridFunctions>
  struct geo_rhs {

    template <class GF>
    using LocalFunction = std::decay_t<decltype(localFunction(Dune::resolveRef(std::declval<GF const&>())))>;

    std::tuple<GridFunctions...> gridFcts_;
    std::tuple<LocalFunction<GridFunctions>...> localFcts_;

    geo_rhs(GridFunctions const&... gridFcts) 
      : gridFcts_(gridFcts...)
      , localFcts_(localFunction(Dune::resolveRef(gridFcts))...)
       {}
  };

} // end namespace tag


template <class... GridFunctions>
class GeoRhsOperator
{
  using Self = GeoRhsOperator;

  template <class GF>
  using LocalFunction = std::decay_t<decltype(localFunction(Dune::resolveRef(std::declval<GF const&>())))>;

public:
  GeoRhsOperator(tag::geo_rhs<GridFunctions...> tag) 
    : gridFcts_(tag.gridFcts_)
    , localFcts_(tag.localFcts_)
  {}


  template <class CG, class Node, class Quad, class LocalFct, class Vec>
  void assemble(CG const& contextGeo, Node const& node, Quad const& /*quad*/,
                LocalFct const& localFct, Vec& elementVector) const
  {
    using namespace Dune::Indices;

    static_assert(CG::dim == 2);
    static_assert(CG::dow == 3);

    auto const& S1Node0 = node.child(_0);
    auto const& S2Node0 = node.child(_1);
    auto const& S3Node0 = node.child(_2);
    auto const& HHNode0 = node.child(_3);

    std::size_t numSLocalFE = S1Node0.child(0).finiteElement().size();
    std::size_t numHHLocalFE = HHNode0.child(0).finiteElement().size();

    using GlobalCoordinate = typename CG::Geometry::GlobalCoordinate;
    using T = typename CG::Geometry::ctype;

    auto const& element = contextGeo.element();
    auto const& geometry = contextGeo.geometry().impl();
    auto const& quad = Dune::QuadratureRules<T, 2>::rule(element.type(), 20);

    for (std::size_t iq = 0; iq < quad.size(); ++iq)
    {
      auto const& qp = quad[iq];

      // The transposed inverse Jacobian of the map from the reference element to the element
      auto const& jacobian = contextGeo.geometry().jacobianInverseTransposed(qp.position());

      // The multiplicative factor in the integral transformation formula
      auto const dS = contextGeo.geometry().integrationElement(qp.position()) * qp.weight();

      // Compute the shape function gradients on the real element
      auto const& SShapeValues = S1Node0.child(0).localBasisValuesAt(qp.position());
      auto const& HHShapeValues = HHNode0.child(0).localBasisValuesAt(qp.position());

      auto const& SShapeGradients = S1Node0.child(0).localBasisJacobiansAt(qp.position());
      auto const& HHShapeGradients = HHNode0.child(0).localBasisJacobiansAt(qp.position());
      std::vector<GlobalCoordinate> SGradients(SShapeGradients.size());
      std::vector<GlobalCoordinate> HHGradients(HHShapeGradients.size());
      for (std::size_t i = 0; i < SGradients.size(); ++i)
        jacobian.mv(SShapeGradients[i][0], SGradients[i]);
      for (std::size_t i = 0; i < HHGradients.size(); ++i)
        jacobian.mv(HHShapeGradients[i][0], HHGradients[i]);

      auto const exprValue = localFct(qp.position());

      // Compute the actual vector entries

      auto const nh = geometry.normal(qp.position());

      // -<n1, div(delta1)>
      for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t i = 0; i < numSLocalFE; ++i) {
        std::size_t local_i = S1Node0.child(k).localIndex(i);
        elementVector[local_i] -= nh[0] * SGradients[i][k]  * dS;
        }
      }

      // -<n2, div(delta2)>
      for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t i = 0; i < numSLocalFE; ++i) {
        std::size_t local_i = S2Node0.child(k).localIndex(i);
        elementVector[local_i] -= nh[1] * SGradients[i][k]  * dS;
        }
      }

      // -<n3, div(delta3)>
      for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t i = 0; i < numSLocalFE; ++i) {
        std::size_t local_i = S3Node0.child(k).localIndex(i);
        elementVector[local_i] -= nh[2] * SGradients[i][k]  * dS;
        }
      }

      // - <grad X, grad gamma>
      for (std::size_t i = 0; i < numHHLocalFE; ++i) {
        for (int r = 0; r < 3; ++r) {
          for (int s = 0; s < 3; ++s) {
            std::size_t local_i = HHNode0.child(s).localIndex(i);
            elementVector[local_i] -= exprValue[s][r] * HHGradients[i][r] * dS;
          }
        }
      }
    }
  }

  private:
  std::tuple<GridFunctions...> gridFcts_;
  mutable std::tuple<LocalFunction<GridFunctions>...> localFcts_;
};

template <class... GridFunctions, class LC>
struct GridFunctionOperatorRegistry<tag::geo_rhs<GridFunctions...>, LC>
{
  static constexpr int degree = 1;
  using type = GeoRhsOperator<GridFunctions...>;
};

/**
 * @}
 **/

} // end namespace AMDiS
