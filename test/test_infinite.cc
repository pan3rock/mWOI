#include "catch.hpp"
#include "mwtrans.hpp"
#include "timer.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <fmt/format.h>
#include <functional>

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

TEST_CASE("infinite integral 1", "[mW transform Integral]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  MWtransInt mwt(0, 0, abserr, referr);
  double ret = mwt.perform(func_1);
  double exact = -1.0;
  REQUIRE(ret == Approx(exact));
}

TEST_CASE("infinite integral 2", "[mW transform Integral]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  MWtransInt mwt(0, 0, abserr, referr);
  double ret = mwt.perform(func_2);
  double exact = 0.421024438240708333;
  REQUIRE(ret == Approx(exact));
}

TEST_CASE("infinite integral 3", "[mW transform Integral]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  MWtransInt mwt(0, 0, abserr, referr);
  double ret = mwt.perform(func_3);
  double exact = 0.421024438240708333;
  REQUIRE(ret == Approx(exact));
}

TEST_CASE("infinite integral 4", "[mW transform Integral]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  MWtransInt mwt(0, 0, abserr, referr);
  double ret = mwt.perform(func_4);
  double exact = 1.0;
  REQUIRE(ret == Approx(exact));
}