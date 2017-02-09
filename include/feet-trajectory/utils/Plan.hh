#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace feettrajectory
{
class Plan
{
 public:
  Plan(int boxBelow, int boxAbove)
  {
    boxBelow_ = boxBelow;
    boxAbove_ = boxAbove;
  }

  const int& boxBelow() const { return boxBelow_; };
  const int& boxAbove() const { return boxAbove_; };

 private:
  int boxBelow_, boxAbove_;
};

} /* feettrajectory */

