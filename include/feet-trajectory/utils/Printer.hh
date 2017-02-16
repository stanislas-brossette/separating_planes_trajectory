#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include <manifolds/Point.h>

#include <feet-trajectory/BoxTrajProblem.hh>
#include <feet-trajectory/BoxesHullTrajProblem.hh>

namespace feettrajectory
{
void print(const std::string& fileName, const BoxTrajProblem& pb,
           const mnf::Point& xStar);
void print(const std::string& fileName, const BoxesHullTrajProblem& pb,
           const mnf::Point& xStar);
void print(const std::string& fileName, const Eigen::Vector3d& bSize,
           const Eigen::Vector3d& oSize, const Eigen::Vector3d& oPos,
           const Eigen::Vector3d& t, const double& d, const Eigen::Vector3d& n);
void printAllIterations(const std::string& fileName,
                        const BoxesHullTrajProblem& pb, const mnf::Point& xStar,
                        const std::string& folder);
std::vector<Eigen::VectorXd> parseX(const std::string& file, const mnf::Point& x);
} /* feettrajectory */

