#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace feettrajectory
{
class PlanForHull
{
  // TODO: plans should be linked to boxes by references or some king of
  // pointers instead of indexes
 public:
  PlanForHull(int boxBelow, int box0Above, int box1Above)
      : boxBelow_(boxBelow), box0Above_(box0Above), box1Above_(box1Above)
  {
  }

  const int& boxBelow() const { return boxBelow_; };
  const int& box0Above() const { return box0Above_; };
  const int& box1Above() const { return box1Above_; };

 private:
  int boxBelow_, box0Above_, box1Above_;
};

} /* feettrajectory */

