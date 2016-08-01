#pragma once

#include <iostream>
#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

#include <manifolds/defs.h>
#include <manifolds/RealSpace.h>
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>
#include <manifolds/CartesianProduct.h>

#include <pgsolver/solver/Problem.h>
#include <pgsolver/utils/defs.h>

#include <cube-stacks/CubeStackProblemOnR.hh>
#include <cube-stacks/functions/TotalCost.hh>
#include <cube-stacks/functions/CubeAboveFixedPlanRoboptim.hh>
#include <cube-stacks/functions/CubeAbovePlanRoboptim.hh>
#include <cube-stacks/functions/Norm1Quaternion.hh>
#include <cube-stacks/functions/Norm1Vector.hh>
#include <cube-stacks/utils/ProblemConfig.hh>

namespace cubestacks
{

class RoboptimCubeStackProblem
{
 public:
  RoboptimCubeStackProblem(const mnf::Manifold& M, const std::string& configPath);
  virtual ~RoboptimCubeStackProblem();

 private:
  ProblemConfig config_;
  CubeStackProblemOnR myProb_;
  std::vector<boost::shared_ptr<CubeAboveFixedPlanRoboptim>> cubeAboveFixedPlanCstrs_; 
  std::vector<boost::shared_ptr<CubeAbovePlanRoboptim>> cubeAbovePlanCstrs_; 
  std::vector<boost::shared_ptr<Norm1Quaternion>> Norm1Quaternions_; 
  std::vector<boost::shared_ptr<Norm1Vector>> Norm1Vectors_; 
  boost::shared_ptr<TotalCost> f_;
};
}
