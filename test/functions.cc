#include "functions.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <complex>
#include <fmt/format.h>

using complex_d = std::complex<double>;
using namespace std::literals::complex_literals;

using boost::math::cyl_bessel_j;
using std::exp;
using std::log;
using std::pow;
using std::sqrt;

double func_1(double x) {
  // \int_0^infty f = -1
  return pow(x, 2) * cyl_bessel_j(0.0, x);
}

double func_2(double x) {
  // 0.421024438240708333
  return 0.5 * log(1.0 + pow(x, 2)) * cyl_bessel_j(1.0, x);
}

double func_3(double x) {
  // 0.421024438240708333
  return x / (1.0 + pow(x, 2)) * cyl_bessel_j(0.0, x);
}

double func_4(double x) {
  // 1.0
  return (1.0 - exp(-x)) / (x * log(1 + sqrt(2))) * cyl_bessel_j(0.0, x);
}

complex_d func_1c(double x) {
  // \int_0^infty f = -1 + 1.0i
  return func_1(x) + 1.0i * func_4(x);
}

complex_d func_simple(double x) { return x + 1.0i * pow(x, 3); }
