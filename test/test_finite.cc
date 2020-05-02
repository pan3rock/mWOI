#include "catch.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <functional>
#include <iostream>

using boost::math::cyl_bessel_j;
using namespace boost::math::quadrature;
using namespace std::placeholders;

double f1(double x, double r) {
  return 1.0 / (pow(x, 2) + 1.0) * cyl_bessel_j(0.0, r * x) * x +
         1.0 / (pow(x, 2) + 2.0) * cyl_bessel_j(1.0, r * x) * x;
}

double f2(double x) { return pow(x, 2); }

TEST_CASE("Guass-Kronrod (f1)", "[mW transform Integral]") {
  double kmin = 1.0;
  double kmax = 2.0;
  double r = 100.0;

  std::function<double(double)> func = std::bind(f1, _1, r);
  double error = 0;
  double estimate =
      gauss_kronrod<double, 61>::integrate(func, kmin, kmax, 0, 0, &error);

  std::cout << "error: " << error << std::endl;
  REQUIRE(estimate == Approx(0.000289257995979269));
}

TEST_CASE("Guass-Kronrod (f2)", "[mW transform Integral]") {
  double kmin = 1.0;
  double kmax = 2.0;

  double error;
  double estimate =
      gauss_kronrod<double, 15>::integrate(f2, kmin, kmax, 0, 0, &error);

  REQUIRE(estimate == Approx(2.3333333333333335));
}