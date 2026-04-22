// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SCALAR_BASIS_HH
#define DUNE_SCALAR_BASIS_HH

#include <cassert>
#include <cstddef>
#include <optional>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/functions/functionspacebases/nodes.hh>


namespace Dune {
namespace Functions {

// forward declarations
template<class GV>
class ScalarNode;

template<class GV>
class ScalarPreBasis;



template<class GV>
class ScalarPreBasis
{
public:

  using Node = ScalarNode<GV>;

  using GridView = GV;
  using size_type = std::size_t;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  ScalarPreBasis(const GridView& gv)
    : gridView_(gv)
  {}

  /// \brief Initialize the global indices. No initialization necessary, since
  /// there is only 1 global index.
  void initializeIndices() {}

  /// \brief Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return gridView_;
  }

  /// \brief Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /// \brief Create tree node
  Node makeNode() const
  {
    return Node{};
  }

  /// \brief The scalar basis represents a 1 dimensional space, so there is
  /// only one index and that is 0
  template<class It>
  It indices(const Node& /*node*/, It it) const
  {
    *it++ = {{ 0u }};
    return it;
  }

  /// \brief Same as size(prefix) with empty prefix
  size_type size() const
  {
    return 1;
  }

  /// \brief Return number of possible values for next position in multi index
  template<class SizePrefix>
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  /// \brief Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return size();
  }

  /// \brief Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return 1;
  }

protected:
  GridView gridView_;
};


template<class GV>
class ScalarNode
    : public LeafBasisNode
{
  using Base = LeafBasisNode;

public:

  // required typedefs
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = P0LocalFiniteElement<typename GV::ctype,double,GV::dimension>;

  /// \brief Return current element, assert if unbound
  const Element& element() const
  {
    assert(element_);
    return *element_;
  }

  /// \brief Return the current bound local finite-element
  const FiniteElement& finiteElement() const
  {
    assert(finiteElement_);
    return *finiteElement_;
  }

  /// \brief Bind to element `e`. Stores a pointer to the element and
  /// creates a local finite-element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_.emplace(element_->type());
    this->setSize(1);
  }

protected:
  std::optional<FiniteElement> finiteElement_;
  const Element* element_ = nullptr;
};


namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a scalar pre-basis.
 * See \ref ScalarBasis for details about a scalar basis.
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
inline auto scalar()
{
  return [](const auto& gridView) {
    return ScalarPreBasis{gridView};
  };
}

} // end namespace BasisFactory


/**
 * \brief A scalar basis represents a 1 dimensional function space `space{1}`
 * of constants. Although bound to a GridView `GV` it is essential independent
 * of any grid elements. On all bound elements the global index of any local DOF
 * is the index 0.
 *
 * The local finite-element is of type \ref P0LocalFiniteElement, i.e. a polynomial
 * basis of order 0.
 *
 * A ScalarBasis can be created with the factory \ref scalar().
 *
 * See the following classes for details:
 * - \ref ScalarPreBasis
 * - \ref ScalarNode
 */
template<class GV>
using ScalarBasis = DefaultGlobalBasis<ScalarPreBasis<GV> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_SCALAR_BASIS_HH
