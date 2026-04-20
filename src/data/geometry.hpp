#pragma once

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace AMDiS {

// Hessian of levelset function phi
struct D2Phi
{
  const double c;
  const double d;

  D2Phi(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  Dune::FieldMatrix<double,3,3> operator()(Dune::FieldVector<double,3> const& X) const;
};

// Gradient of levelset function phi
struct DPhi
{
  const double c;
  const double d;

  DPhi(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const;

  friend D2Phi derivative(DPhi dphi)
  {
    return {dphi.c,dphi.d};
  }
};

// Levelset function phi
struct Phi
{
  const double c;
  const double d;

  Phi(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  double operator()(Dune::FieldVector<double,3> const& X) const;

  friend DPhi derivative(Phi phi)
  {
    return {phi.c,phi.d};
  }
};

// Normal vector of surface
struct N
{
  const double c;
  const double d;

  N(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const;
};

// Mean curvature tr(H) of surface
struct MeanCurvature
{
  const double c;
  const double d;

  MeanCurvature(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  double operator()(Dune::FieldVector<double,3> const& X) const;
};

// Weingarten map H of surface
struct Weingarten
{
  const double c;
  const double d;

  Weingarten(double c = 0.95, double d = 0.96)
    : c(c)
    , d(d)
  {}

  Dune::FieldMatrix<double,3,3> operator()(Dune::FieldVector<double,3> const& X) const;
};

} // end namespace AMDiS
