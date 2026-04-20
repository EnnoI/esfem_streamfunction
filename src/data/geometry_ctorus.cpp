#include <cmath>
#include <iostream>

#include "geometry.hpp"

namespace AMDiS {

// This data implements the Closed-Torus surface

double Phi::operator()(Dune::FieldVector<double,3> const& X) const
{
  double x = X[0], y = X[1], z = X[2];
  return Dune::power(d*d + x*x + y*y + z*z, 3) - 8*d*d*(y*y + z*z) - c*c*c*c;
}


Dune::FieldVector<double,3> DPhi::operator()(Dune::FieldVector<double,3> const& X) const
{
  using std::sqrt;
  double x = X[0], y = X[1], z = X[2];
  double d_2 = d*d; // 5
  double x_2 = x*x; // 3
  double y_2 = y*y; // 3
  double z_2 = z*z; // 3
  double a4 = d_2+x_2+y_2+z_2; // 3
  double a4_2 = a4*a4; // 3
  return {6*x*a4_2,-16*d_2*y+6*y*a4_2,-16*d_2*z+6*z*a4_2};
}

Dune::FieldMatrix<double,3,3> D2Phi::operator()(Dune::FieldVector<double,3> const& X) const
{
  using std::sqrt;
  double x = X[0], y = X[1], z = X[2];
  double x_2 = x*x; // 13
  double d_2 = d*d; // 14
  double y_2 = y*y; // 13
  double z_2 = z*z; // 13
  double a4 = d_2+x_2+y_2+z_2; // 12
  double a4_2 = a4*a4; // 3
  return {{24*x_2*a4+6*a4_2,24*x*y*a4,24*x*z*a4},{24*x*y*a4,-16*d_2+24*y_2*a4+6*a4_2,24*y*z*a4},{24*x*z*a4,24*y*z*a4,-16*d_2+24*z_2*a4+6*a4_2}};
}

Dune::FieldVector<double,3> N::operator()(Dune::FieldVector<double,3> const& X) const
{
  using std::sqrt;
  double x = X[0], y = X[1], z = X[2];
  double d_2 = d*d; // 20
  double x_2 = x*x; // 15
  double y_2 = y*y; // 12
  double z_2 = z*z; // 12
  double a4 = d_2 + x_2 + y_2 + z_2; // 12
  double a4_2 = a4*a4; // 9
  double a4_4 = a4*a4*a4*a4; // 3
  double a5 = 16*d_2*y - 6*y*a4_2;
  double a5_2 = a5*a5;
  double a6 = 16*d_2*z - 6*z*a4_2;
  double a6_2 = a6*a6;
  return {(6*x*a4_2)/sqrt(36*x_2*a4_4 + a5_2 + a6_2), (-16*d_2*y + 6*y*a4_2)/sqrt(36*x_2*a4_4 + a5_2 + a6_2), (-16*d_2*z + 6*z*a4_2)/sqrt(36*x_2*a4_4 + a5_2 + a6_2)};
}


double MeanCurvature::operator()(Dune::FieldVector<double,3> const& X) const
{
  using std::sqrt;
  double x = X[0], y = X[1], z = X[2];
  double d_2 = d*d; // 886
  double x_2 = x*x; // 682
  double y_2 = y*y; // 576
  double z_2 = z*z; // 576
  double x_3 = x*x*x; // 9
  double a5 = d_2+x_2+y_2+z_2; // 564
  double a5_2 = a5*a5; // 350
  double a5_4 = a5*a5*a5*a5; // 106
  double a5_3 = a5*a5*a5; // 27
  double a9 = 16*d_2*y-6*y*a5_2; // 120
  double a10 = 16*d_2*z-6*z*a5_2; // 120
  double a11 = -16*d_2*y+6*y*a5_2; // 29
  double a12 = -16*d_2*z+6*z*a5_2; // 29
  double a13 = 16*d_2-24*y_2*a5-6*a5_2; // 9
  double a14 = 16*d_2-24*z_2*a5-6*a5_2; // 9
  double a15 = 24*x*y*a5; // 6
  double a16 = 24*x*z*a5; // 6
  double a17 = 24*y*z*a5; // 6
  double a18 = 36*x_2*a5_4; // 4
  double a19 = 24*x_2*a5; // 3
  double a20 = 6*a5_2; // 3
  double a21 = -16*d_2+24*y_2*a5+6*a5_2; // 3
  double a22 = -16*d_2+24*z_2*a5+6*a5_2; // 3
  double a11_2 = a11*a11; // 4
  double a9_2 = a9*a9; // 93
  double a10_2 = a10*a10; // 93
  double a12_2 = a12*a12; // 4
  double a27 = 36*x_2*a5_4+a9_2+a10_2; // 93
  double a28 = 288*x_3*a5_3+72*x*a5_4-48*x*y*a5*a9-48*x*z*a5*a10; // 9
  double a29 = 288*x_2*y*a5_3+2*a13*a9-48*y*z*a5*a10; // 9
  double a30 = 288*x_2*z*a5_3-48*y*z*a5*a9+2*a14*a10; // 9
  double a27_3_2 = a27*sqrt(a27); // 27
  double a27_1_2 = sqrt(a27); // 30
  double a33 = 2*a27_3_2; // 18
  double a34 = 1-a11_2/a27; // 4
  double a35 = 1-a12_2/a27; // 4
  double a36 = 1-a18/a27; // 4
  double a37 = a11*a28; // 3
  double a38 = a12*a28; // 3
  double a39 = -3*x*a5_2*a28; // 3
  double a40 = -3*x*a5_2*a29; // 3
  double a41 = a12*a29; // 3
  double a42 = a11*a29; // 3
  double a43 = -3*x*a5_2*a30; // 3
  double a44 = a11*a30; // 3
  double a45 = a12*a30; // 3
  double a46 = -a37/a33+a15/a27_1_2; // 3
  double a47 = -a38/a33+a16/a27_1_2; // 3
  double a48 = a39/a27_3_2+a19/a27_1_2+a20/a27_1_2; // 3
  double a49 = a40/a27_3_2+a15/a27_1_2; // 3
  double a50 = -a41/a33+a17/a27_1_2; // 3
  double a51 = -a42/a33+a21/a27_1_2; // 3
  double a52 = a43/a27_3_2+a16/a27_1_2; // 3
  double a53 = -a44/a33+a17/a27_1_2; // 3
  double a54 = -a45/a33+a22/a27_1_2; // 3
  return (-6*x*a5_2*a11*(a34*a46-(a11*a12*a47)/a27-(6*x*a5_2*a11*a48)/a27))/a27-(6*x*a5_2*a12*(-((a11*a12*a46)/a27)+a35*a47-(6*x*a5_2*a12*a48)/a27))/a27+a36*((-6*x*a5_2*a11*a46)/a27-(6*x*a5_2*a12*a47)/a27+a36*a48)-(6*x*a5_2*a11*(a36*a49-(6*x*a5_2*a12*a50)/a27-(6*x*a5_2*a11*a51)/a27))/a27-(a11*a12*((-6*x*a5_2*a12*a49)/a27+a35*a50-(a11*a12*a51)/a27))/a27+a34*((-6*x*a5_2*a11*a49)/a27-(a11*a12*a50)/a27+a34*a51)-(6*x*a5_2*a12*(a36*a52-(6*x*a5_2*a11*a53)/a27-(6*x*a5_2*a12*a54)/a27))/a27-(a11*a12*((-6*x*a5_2*a11*a52)/a27+a34*a53-(a11*a12*a54)/a27))/a27+a35*((-6*x*a5_2*a12*a52)/a27-(a11*a12*a53)/a27+a35*a54);
}


Dune::FieldMatrix<double,3,3> Weingarten::operator()(Dune::FieldVector<double,3> const& X) const
{
  using std::sqrt;
  double x = X[0], y = X[1], z = X[2];
  double x_2 = x*x; // 2046
  double d_2 = d*d; // 2658
  double y_2 = y*y; // 1728
  double z_2 = z*z; // 1728
  double x_3 = x*x*x; // 27
  double a5 = d_2+x_2+y_2+z_2; // 1692
  double a5_4 = a5*a5*a5*a5; // 318
  double a5_2 = a5*a5; // 1050
  double a5_3 = a5*a5*a5; // 81
  double a9 = 16*d_2*y-6*y*a5_2; // 360
  double a10 = 16*d_2*z-6*z*a5_2; // 360
  double a11 = -16*d_2*y+6*y*a5_2; // 87
  double a12 = -16*d_2*z+6*z*a5_2; // 87
  double a13 = 16*d_2-24*y_2*a5-6*a5_2; // 27
  double a14 = 16*d_2-24*z_2*a5-6*a5_2; // 27
  double a15 = 24*x*y*a5; // 18
  double a16 = 24*x*z*a5; // 18
  double a17 = 24*y*z*a5; // 18
  double a18 = 36*x_2*a5_4; // 12
  double a19 = 24*x_2*a5; // 9
  double a20 = 6*a5_2; // 9
  double a21 = -16*d_2+24*y_2*a5+6*a5_2; // 9
  double a22 = -16*d_2+24*z_2*a5+6*a5_2; // 9
  double a9_2 = a9*a9; // 279
  double a10_2 = a10*a10; // 279
  double a11_2 = a11*a11; // 12
  double a12_2 = a12*a12; // 12
  double a27 = 36*x_2*a5_4+a9_2+a10_2; // 279
  double a28 = 288*x_3*a5_3+72*x*a5_4-48*x*y*a5*a9-48*x*z*a5*a10; // 27
  double a29 = 288*x_2*y*a5_3+2*a13*a9-48*y*z*a5*a10; // 27
  double a30 = 288*x_2*z*a5_3-48*y*z*a5*a9+2*a14*a10; // 27
  double a27_3_2 = a27*sqrt(a27); // 81
  double a27_1_2 = sqrt(a27); // 90
  double a33 = 2*a27_3_2; // 54
  double a34 = 1-a18/a27; // 12
  double a35 = 1-a11_2/a27; // 12
  double a36 = 1-a12_2/a27; // 12
  double a37 = a11*a28; // 9
  double a38 = a12*a28; // 9
  double a39 = -3*x*a5_2*a28; // 9
  double a40 = -3*x*a5_2*a29; // 9
  double a41 = a12*a29; // 9
  double a42 = a11*a29; // 9
  double a43 = -3*x*a5_2*a30; // 9
  double a44 = a11*a30; // 9
  double a45 = a12*a30; // 9
  double a46 = -a37/a33+a15/a27_1_2; // 9
  double a47 = -a38/a33+a16/a27_1_2; // 9
  double a48 = a39/a27_3_2+a19/a27_1_2+a20/a27_1_2; // 9
  double a49 = a40/a27_3_2+a15/a27_1_2; // 9
  double a50 = -a41/a33+a17/a27_1_2; // 9
  double a51 = -a42/a33+a21/a27_1_2; // 9
  double a52 = a43/a27_3_2+a16/a27_1_2; // 9
  double a53 = -a44/a33+a17/a27_1_2; // 9
  double a54 = -a45/a33+a22/a27_1_2; // 9
  double a55 = -6*x*a5_2*a11*a46; // 3
  double a56 = 6*x*a5_2*a12*a47; // 3
  double a57 = 6*x*a5_2*a12*a50; // 3
  double a58 = 6*x*a5_2*a11*a51; // 3
  double a59 = 6*x*a5_2*a11*a53; // 3
  double a60 = 6*x*a5_2*a12*a54; // 3
  double a61 = a11*a12*a47; // 3
  double a62 = 6*x*a5_2*a11*a48; // 3
  double a63 = -6*x*a5_2*a11*a49; // 3
  double a64 = a11*a12*a50; // 3
  double a65 = -6*x*a5_2*a11*a52; // 3
  double a66 = a11*a12*a54; // 3
  double a67 = a11*a12*a46; // 3
  double a68 = 6*x*a5_2*a12*a48; // 3
  double a69 = -6*x*a5_2*a12*a49; // 3
  double a70 = a11*a12*a51; // 3
  double a71 = -6*x*a5_2*a12*a52; // 3
  double a72 = a11*a12*a53; // 3
  double a73 = a55/a27-a56/a27+a34*a48; // 3
  double a74 = a34*a49-a57/a27-a58/a27; // 3
  double a75 = a34*a52-a59/a27-a60/a27; // 3
  double a76 = a35*a46-a61/a27-a62/a27; // 3
  double a77 = a63/a27-a64/a27+a35*a51; // 3
  double a78 = a65/a27+a35*a53-a66/a27; // 3
  double a79 = a67/a27; // 3
  double a80 = a69/a27+a36*a50-a70/a27; // 3
  double a81 = a71/a27-a72/a27+a36*a54; // 3
  double a82 = -a79+a36*a47-a68/a27; // 3
  return {{a34*a73-(6*x*a5_2*a11*a74)/a27-(6*x*a5_2*a12*a75)/a27,(-6*x*a5_2*a11*a73)/a27+a35*a74-(a11*a12*a75)/a27,(-6*x*a5_2*a12*a73)/a27-(a11*a12*a74)/a27+a36*a75},{a34*a76-(6*x*a5_2*a11*a77)/a27-(6*x*a5_2*a12*a78)/a27,(-6*x*a5_2*a11*a76)/a27+a35*a77-(a11*a12*a78)/a27,(-6*x*a5_2*a12*a76)/a27-(a11*a12*a77)/a27+a36*a78},{a34*a82-(6*x*a5_2*a11*a80)/a27-(6*x*a5_2*a12*a81)/a27,(-6*x*a5_2*a11*a82)/a27+a35*a80-(a11*a12*a81)/a27,(-6*x*a5_2*a12*a82)/a27-(a11*a12*a80)/a27+a36*a81}};
}

} // end namespace AMDiS
