#include "catch.hpp"
#include "mwtrans.hpp"

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

ArrayXd fun1(const ArrayXd &x) { return 1.0 / (pow(x, 2) + 1.0); }

ArrayXd fun2(const ArrayXd &x) { return 1.0 / (pow(x, 2) + 2.0); }

TEST_CASE("Finite integral", "[mW transform Integral]") {
  double kmin = 1.0;
  double kmax = 2.0;
  int nk = 50;
  double r = 100.0;
  auto k = ArrayXd::LinSpaced(nk, kmin, kmax);
  auto f1 = fun1(k);
  auto f2 = fun2(k);

  double estimate = 0;
  std::cout << estimate << std::endl;
  REQUIRE(estimate == Approx(0.000289257995979269));
}