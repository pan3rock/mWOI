#ifndef MWTRANS_H_
#define MWTRANS_H_

#include <Eigen/Dense>
#include <functional>

class MWtransInt {
public:
  MWtransInt(double lb, double v, double abserr, double referr, int max_depth);

  double perform(std::function<double(double)> functor);

private:
  const int nzero_ = 100;
  const int niter_ = 100;
  const unsigned max_depth_;

  Eigen::ArrayXd zeros_;
  double lb_, ub_, v_;
  const double abserr_;
  const double referr_;
};

#endif