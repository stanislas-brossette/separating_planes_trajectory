#pragma once
#include <iostream>

namespace cubestacks
{
class ProblemInfo
{
 public:
  ProblemInfo(int nCubes) : nCubes_(nCubes)
  {
    nPlanes_ = nCubes_ * (nCubes_ - 1) / 2;
    totalDim_ = 7*nCubes_ + 4*nPlanes_;
  }
  int nCubes_;
  int nPlanes_;
  int totalDim_;
};
} /* cubestacks */ 
