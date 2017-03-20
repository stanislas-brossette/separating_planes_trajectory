#pragma once

#include <iostream>
#include <Eigen/Core>

namespace feettrajectory
{
class QP
{
 public:
  QP();
  QP(const long& n, const long& m);
  void setDimensions(const long& n, const long& m);
  virtual ~QP();

 protected:
  long n_;             // Size variables
  long m_;             // Size constraints
  Eigen::MatrixXd A_;  // Quadratic cost
  Eigen::MatrixXd C_;  // Constraint
  Eigen::VectorXd c_;  // linear cost
  Eigen::VectorXd l_;  // lower bound
  Eigen::VectorXd u_;  // upper bound
};
} /* feettrajectory */
