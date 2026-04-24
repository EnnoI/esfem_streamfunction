#include <src/data/geometry.hpp>
#include <dune/curvedgrid/geometries/implicitsurface.hh>


struct Surface
{
  Surface(double c = 0.95, double d = 0.96)
    : phi{c,d}
    , projection{phi, 100}
    , normal{c,d}
    , mean_curvature{c,d}
    , weingarten{c,d}
  {}

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const
  {
    return projection(X);
  }

  AMDiS::Phi phi;
  Dune::ImplicitSurfaceProjection<AMDiS::Phi> projection;

  AMDiS::N normal;
  AMDiS::MeanCurvature mean_curvature;
  AMDiS::Weingarten weingarten;
};

struct Sphere
{

  Sphere()
    : normal{}
    , mean_curvature{}
    , weingarten{}
    , projection{}
  {}

  struct Projection
  {
    Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const
    {
      return X / X.two_norm();
    }
  };

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const
  {
    return X / X.two_norm();
  };

  struct Normal
  {
    Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& X) const
    {
      return X / X.two_norm();
    }
  };

  struct MeanCurvature
  {
    double operator()(Dune::FieldVector<double,3> const& X) const
    {
      return -2.0; // curvature of unit sphere
    }
  };

  struct Weingarten
  {
    Dune::FieldMatrix<double,3,3> operator()(Dune::FieldVector<double,3> const& X) const
    {
      return 1; // curvature of unit sphere
    }
  };

  Normal normal;
  MeanCurvature mean_curvature;
  Weingarten weingarten;
  Projection projection;

};