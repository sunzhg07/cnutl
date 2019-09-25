
#include "coupling_coeff.h"

using std::abs;
using std::cout;
using std::endl;
using std::max;
using std::min;

double Coupling::add(double a, double b) {
  double val = a + b;
  return val;
}
double Coupling::SixJ(int j2a, int j2b, int J2ab, int j2c, int J2, int J2bc) {
  double val = gsl_sf_coupling_6j(j2a, j2b, J2ab, j2c, J2, J2bc);
  return val;
}
double Coupling::NinJ(int j2a, int j2b, int J2ab, int j2c, int j2d, int J2cd,
                      int J2ac, int J2bd, int J2) {
  double val =
      gsl_sf_coupling_9j(j2a, j2b, J2ab, j2c, j2d, J2cd, J2ac, J2bd, J2);

  // cout << j2a << " " << j2b << " " << J2ab << " " << j2c << " " << j2d << " "
  //      << J2cd << " " << J2ac << " " << J2bd << " " << J2 << endl;
  // cout << " Val_cpp = " << val << endl;
  return val;
}
double Coupling::factorial(double n) {
  if (n < 0)
    return -100.;

  return (n == 1. || n == 0.) ? 1. : factorial(n - 1) * n;
}


double Coupling::CGcoeff(int j2a, int j2b, int J2, int m2a, int m2b, int M2) {
  double result_3j, cg;
  int M2_x = M2 * (-1);
  int phase;
  if ((m2a + m2b) != M2)
    return 0.0;
  if (std::abs(j2a - j2b) > J2 || (j2a + j2b) < J2)
    return 0.0;
  result_3j = gsl_sf_coupling_3j(j2a, j2b, J2, m2a, m2b, M2_x);
  phase = 1 - 2 * (int(std::abs((j2a - j2b + M2) / 2)) % 2);
  cg = result_3j * std::sqrt(J2 + 1) * phase;
  return cg;
}
