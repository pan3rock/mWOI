#include "catch.hpp"
#include "functions.hpp"
#include "mwtrans.hpp"
#include "timer.hpp"

#include <fmt/format.h>
#include <functional>

TEST_CASE("Case 1", "[test_infinite]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  unsigned max_depth = 5;
  MWtransInt mwt(0, 0, abserr, referr, max_depth);
  double ret = mwt.perform(func_1);
  double exact = -1.0;
  REQUIRE(ret == Approx(exact));
}

TEST_CASE("Case 2", "[test_infinite]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  unsigned max_depth = 5;
  MWtransInt mwt(0, 0, abserr, referr, max_depth);
  double ret = mwt.perform(func_2);
  double exact = 0.421024438240708333;
  REQUIRE(ret == Approx(exact));
}

TEST_CASE("Case 3", "[test_infinite]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  unsigned max_depth = 5;
  MWtransInt mwt(0, 0, abserr, referr, max_depth);
  double ret = mwt.perform(func_3);
  double exact = 0.421024438240708333;
  REQUIRE(ret == Approx(exact));
}

TEST_CASE("Case 4", "[test_infinite]") {
  double abserr = 1.0e-9;
  double referr = 1.0e-7;
  unsigned max_depth = 5;
  MWtransInt mwt(0, 0, abserr, referr, max_depth);
  double ret = mwt.perform(func_4);
  double exact = 1.0;
  REQUIRE(ret == Approx(exact));
}

TEST_CASE("benchmark abs", "[test_infinite]") {
  BENCHMARK("abs 1.0e-9") {
    double abserr = 1.0e-9;
    double referr = 1.0e-6;
    unsigned max_depth = 5;
    MWtransInt mwt(0, 0, abserr, referr, max_depth);
    double ret = mwt.perform(func_4);
    double exact = 1.0;
    REQUIRE(ret == Approx(exact));
  };
  BENCHMARK("abs 1.0e-7") {
    double abserr = 1.0e-7;
    double referr = 1.0e-6;
    unsigned max_depth = 5;
    MWtransInt mwt(0, 0, abserr, referr, max_depth);
    double ret = mwt.perform(func_4);
    double exact = 1.0;
    REQUIRE(ret == Approx(exact));
  };
  BENCHMARK("abs 1.0e-5") {
    double abserr = 1.0e-5;
    double referr = 1.0e-6;
    unsigned max_depth = 5;
    MWtransInt mwt(0, 0, abserr, referr, max_depth);
    double ret = mwt.perform(func_4);
    double exact = 1.0;
    REQUIRE(ret == Approx(exact));
  };
}

TEST_CASE("benchmark max_depth", "[test_infinite]") {
  BENCHMARK("max_depth 5") {
    double abserr = 1.0e-9;
    double referr = 1.0e-6;
    unsigned max_depth = 5;
    MWtransInt mwt(0, 0, abserr, referr, max_depth);
    double ret = mwt.perform(func_4);
    double exact = 1.0;
    REQUIRE(ret == Approx(exact));
  };
  BENCHMARK("max_depth 10") {
    double abserr = 1.0e-9;
    double referr = 1.0e-6;
    unsigned max_depth = 10;
    MWtransInt mwt(0, 0, abserr, referr, max_depth);
    double ret = mwt.perform(func_4);
    double exact = 1.0;
    REQUIRE(ret == Approx(exact));
  };
  BENCHMARK("max_depth 15") {
    double abserr = 1.0e-9;
    double referr = 1.0e-6;
    unsigned max_depth = 15;
    MWtransInt mwt(0, 0, abserr, referr, max_depth);
    double ret = mwt.perform(func_4);
    double exact = 1.0;
    REQUIRE(ret == Approx(exact));
  };
}