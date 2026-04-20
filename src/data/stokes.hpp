#pragma once

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace AMDiS {

// Velocity u
struct U
{
  const double c;
  const double d;

  U(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const;
};

// Pressure p
struct P
{
  const double c;
  const double d;

  P(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  double operator()(Dune::FieldVector<double,3> const& X) const;
};

// rhs-function f
struct F
{
  const double c;
  const double d;
  std::vector<double> parameters;

  F(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {
    initParameters();
  }

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const;

  void test() const;
  void initParameters();
};

// vorticity psi
struct Psi
{
  const double c;
  const double d;

  Psi(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  double operator()(Dune::FieldVector<double,3> const& X) const;
};

// laplace-beltrami of vorticity psi
struct LaplacePsi
{
  const double c;
  const double d;

  LaplacePsi(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  double operator()(Dune::FieldVector<double,3> const& X) const;
};

} // end namespace AMDiS
