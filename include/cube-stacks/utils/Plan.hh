#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace cubestacks
{

class Plan
{
 public:
  Plan(size_t cubeBelow, size_t cubeAbove)
  {
    cubeBelow_ = cubeBelow;
    cubeAbove_ = cubeAbove;
  }

  const size_t& cubeBelow() const { return cubeBelow_; };
  const size_t& cubeAbove() const { return cubeAbove_; };

 private:
  size_t cubeBelow_, cubeAbove_;
};
  
} /* cubestacks */ 

