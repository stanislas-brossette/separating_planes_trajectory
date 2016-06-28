#pragma once

#include <iostream>
#include <Eigen/Core>

#include <manifolds/defs.h>
#include <pgsolver/solver/Problem.h>
#include <manifolds/RealSpace.h>
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>
#include <manifolds/CartesianProduct.h>

#include <pgsolver/utils/defs.h>

#include <cube-stacks/utils/Cube.hh>
#include <cube-stacks/utils/Plan.hh>
#include <cube-stacks/utils/ProblemConfig.hh>
#include <cube-stacks/functions/CubeAbovePlan.hh>
#include <cube-stacks/functions/CubeAboveFixedPlan.hh>

namespace cubestacks
{

class CubeStackProblemOnManifold : public pgs::Problem
{
 public:
  CubeStackProblemOnManifold(const mnf::Manifold& M, const std::string& configPath);

  static mnf::CartesianProduct* buildManifold(Index nCubes);
  Eigen::VectorXd findInitPoint();

  void getTangentLB(RefVec out) const;
  void getTangentUB(RefVec out) const;

  void evalObj(double& out) const;
  void evalObjDiff(RefMat out) const;

  void evalLinCstr(RefVec out, size_t i) const;
  void evalLinCstrDiff(RefMat out, size_t i) const;
  void getLinCstrLB(RefVec out, size_t i) const;
  void getLinCstrUB(RefVec out, size_t i) const;

  void evalNonLinCstr(RefVec out, size_t i) const;
  void evalNonLinCstrDiff(RefMat out, size_t i) const;
  void getNonLinCstrLB(RefVec out, size_t i) const;
  void getNonLinCstrUB(RefVec out, size_t i) const;

  size_t numberOfCstr() const;
  Index linCstrDim(size_t i) const;
  Index nonLinCstrDim(size_t i) const;

  void fileForMatlab(std::string fileName, const mnf::Point& x) const;

  std::string getCstrName(const size_t i) const;

 public:
  size_t nCubes_;
  size_t nPlans_;

 private:
  std::vector<Cube> cubes_;
  std::vector<Plan> plans_;
  std::vector<CubeAboveFixedPlan> cubeAboveFixedPlanCstrs_;
  std::vector<CubeAbovePlan> cubeAbovePlanFcts_;
  std::vector<std::string> cstrNames_;
  double distZPlus_, distXPlus_, distXMinus_, distYPlus_, distYMinus_;
  Eigen::Vector3d normalZPlus_, normalXPlus_, normalXMinus_, normalYPlus_,
      normalYMinus_;
  ProblemConfig config_;

  // buffers
  mutable Eigen::MatrixXd outRepObjDiff_;
  //mutable Eigen::MatrixXd outRepLinCstrDiff_;
  mutable Eigen::MatrixXd outRep_;
};
}
