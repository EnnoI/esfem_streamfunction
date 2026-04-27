#include <src/data/geometry.hpp>
#include <dune/curvedgrid/geometries/implicitsurface.hh>
#include <boost/math/special_functions/spherical_harmonic.hpp>


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


struct SphericalHarmonic {

   SphericalHarmonic(unsigned const l, int const m, double const factor, double const radius)
    : mean_curvature{}
    , l_(l)
    , m_(m)
    , factor_(factor)
    , radius_(radius)
     {}

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& x) const
  {
    double theta{std::acos(x[2]/x.two_norm())};
    double phi{std::abs(x[2]-x.two_norm())/x.two_norm() > 1.e-4 ? std::atan2(x[1],x[0]) : 0.};

    double perturb{0.};
    if (m_ < 0) {
      perturb = radius_ + factor_ * std::sqrt(2.) * std::pow(-1,m_) * boost::math::spherical_harmonic_i(l_, -m_, theta, phi);
    } else if (m_ == 0) {
      perturb = radius_ + factor_ * boost::math::spherical_harmonic_r(l_, m_, theta, phi);
    } else {
      perturb = radius_ + factor_ * std::sqrt(2.) * std::pow(-1,m_) * boost::math::spherical_harmonic_r(l_, m_, theta, phi);
    }

    return Dune::FieldVector<double, 3>{perturb * x[0]/x.two_norm(),
                                        perturb * x[1]/x.two_norm(),
                                        perturb * x[2]/x.two_norm()};
  };

  struct MeanCurvature
  {
    double operator()(Dune::FieldVector<double,3> const& X) const
    {
      return -2.0; // curvature of unit sphere
    }
  };

  MeanCurvature mean_curvature;

  private:
  unsigned const l_;
  int const m_;
  double const factor_;
  double const radius_;
};

struct PerturbedSphere {

  PerturbedSphere()
    : mean_curvature{}
     {}

  Dune::FieldVector<double,3> operator()(Dune::FieldVector<double,3> const& x) const
  {

    double theta{std::acos(x[2]/x.two_norm())};
    double phi{std::abs(x[2]-x.two_norm())/x.two_norm() > 1.e-4 ? std::atan2(x[1],x[0]) : 0.};

    double constexpr r0 = 0.4;
    double R = 1.0 + r0 * std::cos(phi) * std::sin(3*theta);
    return R * x / x.two_norm();
  };

  struct MeanCurvature
  {
    double operator()(Dune::FieldVector<double,3> const& X) const
    {
      return -2.0; // curvature of unit sphere
    }
  };

  MeanCurvature mean_curvature;
};