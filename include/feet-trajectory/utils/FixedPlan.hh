#pragma once

namespace feettrajectory
{
class FixedPlan
{
 public:
  FixedPlan(double d, Eigen::Vector3d n) : d_(d), normal_(n) {}
  virtual ~FixedPlan() {}
  const double& d() const { return d_; }
  const Eigen::Vector3d& normal() const { return normal_; }

 private:
  double d_;
  Eigen::Vector3d normal_;
};
} /* feettrajectory */
