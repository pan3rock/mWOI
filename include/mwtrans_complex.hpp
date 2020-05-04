#ifndef MWTRANS_H_
#define MWTRANS_H_

#include <Eigen/Dense>
#include <complex>
#include <functional>

class MWtransIntComplex {
public:
  MWtransIntComplex(double lb, double v, double abserr, double referr);

  std::complex<double>
  perform(std::function<std::complex<double>(double)> functor);

private:
  const int nzero_ = 100;
  const int niter_ = 100;
  const unsigned max_depth_ = 15;

  Eigen::ArrayXd zeros_;
  double lb_, ub_, v_;
  const double abserr_;
  const double referr_;
};

#endif