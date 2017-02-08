#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace feettrajectory
{
class Plan
{
 public:
  Plan(size_t boxBelow, size_t boxAbove)
  {
    boxBelow_ = boxBelow;
    boxAbove_ = boxAbove;
  }

  const size_t& boxBelow() const { return boxBelow_; };
  const size_t& boxAbove() const { return boxAbove_; };

 private:
  size_t boxBelow_, boxAbove_;
};

} /* feettrajectory */

