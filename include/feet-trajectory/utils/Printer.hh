#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include <manifolds/Point.h>

#include <feet-trajectory/BoxTrajProblem.hh>

namespace feettrajectory
{
void print(const std::string& fileName, const BoxTrajProblemOnManifold& pb,
           const mnf::Point& xStar);
void print(const std::string& fileName, const Eigen::Vector3d& bSize,
           const Eigen::Vector3d& oSize, const Eigen::Vector3d& oPos,
           const Eigen::Vector3d& t, const double& d, const Eigen::Vector3d& n);
} /* feettrajectory */

