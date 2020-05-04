#include "catch.hpp"
#include "functions.hpp"
#include "mwtrans_complex.hpp"
#include "timer.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <complex>
#include <fmt/format.h>

using complex_d = std::complex<double>;
using namespace std::literals::complex_literals;
using namespace boost::math::quadrature;

bool is_equal(complex_d c1, complex_d c2, double eps = 1.0e-8) {
  return std::abs(c1 - c2) < eps;
}

TEST_CASE("Case simple", "[test_complex]") {
  double error;
  double a{0};
  double b{1};
  unsigned int max_depth = 0;
  double tolerance = 0;
  complex_d ret = gauss_kronrod<double, 61>::integrate(
      func_simple, a, b, max_depth, tolerance, &error);
  complex_d exact{0.5, 0.25};
  REQUIRE(is_equal(ret, exact));
}

TEST_CASE("Case 1c", "[test_complex]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  MWtransIntComplex mwt(0, 0, abserr, referr);
  complex_d ret = mwt.perform(func_1c);
  complex_d exact = -1.0 + 1.0i;
  std::cout << ret << std::endl;
  REQUIRE(is_equal(ret, exact));
}
