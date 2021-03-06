#pragma once

#include <iostream>
#include <Eigen/Core>

#include <manifolds/defs.h>
#include <manifolds/RealSpace.h>
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>
#include <manifolds/CartesianProduct.h>

#include <pgsolver/solver/Problem.h>
#include <pgsolver/utils/defs.h>

#include <feet-trajectory/utils/Box.hh>
#include <feet-trajectory/utils/Plan.hh>
#include <feet-trajectory/utils/ProblemConfig.hh>
#include <feet-trajectory/functions/BoxAbovePlan.hh>

namespace feettrajectory
{
class BoxTrajProblem: public pgs::Problem
{
 public:
  BoxTrajProblem(const mnf::Manifold& M,
                           const std::string& configPath);
  virtual ~BoxTrajProblem();

  static mnf::CartesianProduct* buildManifold(const Index& nBoxes,
                                              const Index& nObstacles);
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

  // void fileForMatlab(std::string fileName, const mnf::Point& x) const;

  std::string getCstrName(const size_t i) const;

  // getters
  const size_t& nBoxes() const { return nBoxes_; }
  const size_t& nPlans() const { return nPlans_; }
  const size_t& nObstacles() const { return nObstacles_; }
  const std::vector<Box>& boxes() const { return boxes_; }
  const std::vector<Box>& obstacles() const { return obstacles_; }
  const std::vector<Plan>& plans() const { return plans_; }
  const Eigen::Vector3d& initPos() const { return initPos_; }
  const Eigen::Vector3d& finalPos() const { return finalPos_; }
  const Eigen::Vector3d& boxSize() const { return boxSize_; }

 private:
  size_t nBoxes_;
  size_t nPlans_;
  size_t nObstacles_;

  std::vector<Box> boxes_;
  std::vector<Box> obstacles_;
  std::vector<Plan> plans_;
  // std::vector<CubeAboveFixedPlan> cubeAboveFixedPlanCstrs_;
  std::vector<BoxAbovePlan> boxAbovePlanFcts_;
  std::vector<BoxAbovePlan> obstacleAbovePlanFcts_;
  std::vector<std::string> cstrNames_;
  ProblemConfig config_;

  Eigen::Vector3d initPos_;
  Eigen::Vector3d finalPos_;
  Eigen::Vector3d boxSize_;

  // buffers
  mutable Eigen::MatrixXd outRepObjDiff_;
  // mutable Eigen::MatrixXd outRepLinCstrDiff_;
  mutable Eigen::MatrixXd outRep_;
};

} /* feettrajectory */
