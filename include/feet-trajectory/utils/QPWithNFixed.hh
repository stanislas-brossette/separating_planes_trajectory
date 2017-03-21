#pragma once

#include <iostream>
#include <feet-trajectory/utils/QP.hh>
#include <feet-trajectory/TrajectoryProblem.hh>

namespace feettrajectory
{
class QPWithNFixed : public QP
{
 public:
  QPWithNFixed(const TrajectoryProblem& pb);
  virtual ~QPWithNFixed();
  void formQP(const TrajectoryProblem& pb, const Eigen::VectorXd& normals);

 private:
  size_t dimBox_;
  size_t dimBoxes_;
  size_t dimVar_;
  size_t dimCstr_;
};

} /* feettrajectory */

