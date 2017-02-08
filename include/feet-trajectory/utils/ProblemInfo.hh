#pragma once
#include <iostream>

namespace feettrajectory
{
class ProblemInfo
{
 public:
  ProblemInfo(int nboxes) : nBoxes_(nBoxes)
  {
    nPlanes_ = nBoxes_ * (nBoxes_ - 1) / 2;
    totalDim_ = 7 * nBoxes_ + 4 * nPlanes_;
  }
  int nBoxes_;
  int nPlanes_;
  int totalDim_;
};
} /* feettrajectory */
