#ifndef MWTRANS_H_
#define MWTRANS_H_

#include <Eigen/Dense>

class MWtransInt {
public:
  MWtransInt(double lb, double ub);

private:
  double lb_, ub_;
};

#endif