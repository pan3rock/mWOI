#include "mwtrans.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <fmt/format.h>
#include <iostream>
#include <limits>

using namespace Eigen;
using namespace boost::math::quadrature;
using boost::math::cyl_bessel_j_zero;
using std::abs;

MWtransInt::MWtransInt(double lb, double v, double abserr, double referr)
    : zeros_(nzero_), lb_(lb), v_(v), abserr_(abserr), referr_(referr) {
  int nz = 0;
  int i = 0;
  while (nz < nzero_) {
    double r = cyl_bessel_j_zero(v_, ++i);
    if (r > lb_) {
      zeros_(nz) = r;
      ++nz;
    }
  }
  ub_ = zeros_(0) - lb_;
}

double MWtransInt::MWtransInt::perform(std::function<double(double)> functor) {
  // int_lb^z0 functor
  double error;
  double int1 = gauss_kronrod<double, 61>::integrate(functor, lb_, zeros_(0), 0,
                                                     0, &error);

  // int_z0^\infty functor
  ArrayXXd aM = ArrayXXd::Zero(niter_, niter_);
  ArrayXXd aN = ArrayXXd::Zero(niter_, niter_);
  ArrayXd aW = ArrayXd::Zero(niter_);
  ArrayXd x = ArrayXd::Zero(niter_);

  aW(0) = 1.0;
  aW(1) = 1.0;

  int i0 = 0;
  double left = zeros_(i0++);
  double right = zeros_(i0++);

  double fxs =
      gauss_kronrod<double, 15>::integrate(functor, left, right, 0, 0, &error);
  left = right;
  right = zeros_(i0++);
  x(0) = right;

  int neval = 1;
  int count = 0;
  double diff = std::numeric_limits<double>::max();
  double errorbound = 0.0;

  while (count < niter_ - 3 && diff > errorbound) {
    double phixs = gauss_kronrod<double, 15>::integrate(functor, left, right, 0,
                                                        0, &error);
    if (phixs == 0 || abs(phixs) < 1.0e-100) {
      diff = 0.0;
      break;
    }
    left = right;
    right = zeros_(i0++);
    x(count + 1) = right;
    ++neval;

    aM(count, 0) = fxs / phixs;
    aN(count, 0) = 1.0 / phixs;
    fxs += phixs;

    if (count == 1) {
      double denom = 1.0 / (1.0 / x(0) - 1.0 / x(1));
      aM(0, 0) = (aM(0, 0) - aM(1, 0)) * denom;
      aN(0, 0) = (aN(0, 0) - aN(1, 0)) * denom;
      denom = 1.0 / (1.0 / x(0) - 1.0 / x(2));
      aM(0, 1) = (aM(0, 0) - aM(1, 0)) * denom;
      aN(0, 1) = (aN(0, 0) - aN(1, 0)) * denom;
    } else if (count > 1) {
      for (int i = 1; i <= count; ++i) {
        double denom = 1.0 / (1.0 / x(count - i) - 1.0 / x(count + 1));
        aM(count - i, i) =
            (aM(count - i, i - 1) - aM(count + 1 - i, i - 1)) * denom;
        aN(count - i, i) =
            (aN(count - i, i - 1) - aN(count + 1 - i, i - 1)) * denom;
      }
    }
    aW(count) = aM(0, count) / aN(0, count);
    errorbound = std::max(abserr_, referr_ * abs(aW(count)));

    if (count > 1) {
      diff = std::max(abs(aW(count) - aW(count - 1)),
                      abs(aW(count) - aW(count - 2)));
    }
    if (std::isnan(diff)) {
      throw "Warning: Divergent partial sums for mW tranformation.";
    }
    count++;
  }
  if ((count == niter_ - 3) && (diff > errorbound)) {
    throw "Warning: Possible non-convergence of mW transformation";
  }

  double result;
  if (count > 0) {
    result = aW(count - 1);
  } else {
    result = aW(count);
    if (count == 0)
      result = 0;
  }

  return result + int1;
}