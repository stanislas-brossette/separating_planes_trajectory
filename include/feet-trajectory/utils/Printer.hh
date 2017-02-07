#pragma once
#include <iostream>
#include <fstream>
#include <string>


#include <manifolds/Point.h>

#include <feet-trajectory/BoxTrajProblem.hh>

namespace feettrajectory
{
  void print(const std::string& fileName, const BoxTrajProblemOnManifold& pb, const mnf::Point& xStar);
} /* feettrajectory */ 

