#pragma once
#include <iostream>
#include <Eigen/Core>
#include <feet-trajectory/utils/Box.hh>
#include <feet-trajectory/utils/defs.hh>

namespace feettrajectory
{
class FixedBoxPosition
{
 public:
  FixedBoxPosition(const Box& box, const Eigen::Vector3d& targetPos);
  virtual ~FixedBoxPosition();
  // This constraint fixes the position of a given box to a target position
  static void fillLinCstr(const Eigen::Vector3d& targetPos, RefVec lb, RefMat C, RefVec ub);

 private:
  Box box_;
  Eigen::Vector3d targetPos_;
};
} /* feettrajectory */
