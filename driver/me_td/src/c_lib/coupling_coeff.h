#ifndef _COUPLING_COEFF_
#define _COUPLING_COEFF_

#include "gsl/gsl_sf_coupling.h"
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

class Coupling {
public:
  Coupling() { a = 0; }
  double add(double a, double b);
  double SixJ(int j2a, int j2b, int J2ab, int j2c, int J2, int J2bc);
  double NinJ(int j2a, int j2b, int J2ab, int j2c, int j2d, int J2cd, int J2ac,
              int J2bd, int J2);
  double factorial(double n);
  double CGcoeff(int j2a, int j2b, int J2, int m2a, int m2b, int M2);
  
  int a;
};

#endif
