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

typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericFunction<T> Function_t;
typedef roboptim::GenericLinearFunction<T> LinearFunction_t;
typedef roboptim::Solver<T> solver_t;
typedef roboptim::OptimizationLogger<solver_t> logger_t;

struct resultRoboptimCubeStackProblem
{
  Eigen::VectorXd xSol;
  int status;
  double obj_star;
};

class RoboptimCubeStackProblem
{
 public:
  RoboptimCubeStackProblem(const mnf::Manifold& M,
                           const std::string& configPath);
  virtual ~RoboptimCubeStackProblem();
  void init(const Eigen::VectorXd& x0);
  resultRoboptimCubeStackProblem solve(const std::string& solverName);
  void name(const std::string& n);

  size_t nCubes_;
  size_t nPlans_;
  int nIterations_;
  int status_;
  double obj_star_;

  CubeStackProblemOnR& pgsProb() { return pgsProb_; };

  void setCFSQPparameterInt(solver_t& solver, const std::string& s);
  void setCFSQPparameterDouble(solver_t& solver, const std::string& s);
  void setCFSQPparameters(solver_t& solver);

  std::vector<boost::shared_ptr<CubeAboveFixedPlanRoboptim>>
  cubeAboveFixedPlanCstrs()
  {
    return cubeAboveFixedPlanCstrs_;
  }
  std::vector<boost::shared_ptr<CubeAbovePlanRoboptim>> cubeAbovePlanCstrs()
  {
    return cubeAbovePlanCstrs_;
  }
  std::vector<boost::shared_ptr<Norm1Quaternion>> Norm1Quaternions()
  {
    return Norm1Quaternions_;
  }
  std::vector<boost::shared_ptr<Norm1Vector>> Norm1Vectors()
  {
    return Norm1Vectors_;
  }

 private:
  ProblemConfig config_;
  CubeStackProblemOnR pgsProb_;
  IndexManager indexManager_;
  std::vector<boost::shared_ptr<CubeAboveFixedPlanRoboptim>> cubeAboveFixedPlanCstrs_;
  std::vector<boost::shared_ptr<CubeAbovePlanRoboptim>> cubeAbovePlanCstrs_;
  std::vector<boost::shared_ptr<Norm1Quaternion>> Norm1Quaternions_;
  std::vector<boost::shared_ptr<Norm1Vector>> Norm1Vectors_;
  boost::shared_ptr<TotalCost> f_;
  std::string name_;
  std::shared_ptr<solver_t::problem_t> problem_;
};
}

