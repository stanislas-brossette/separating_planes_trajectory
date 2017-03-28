#pragma once

#include <iostream>
#include <Eigen/Core>

#include <feet-trajectory/utils/Box.hh>
#include <feet-trajectory/utils/FixedPlan.hh>
#include <feet-trajectory/utils/PlanForHull.hh>
#include <feet-trajectory/utils/ProblemConfig.hh>
#include <feet-trajectory/functions/BoxAbovePlan.hh>
#include <feet-trajectory/functions/BoxAboveFixedPlan.hh>
#include <feet-trajectory/functions/CostDistance.hh>

namespace feettrajectory
{
class TrajectoryProblem
{
 public:
  TrajectoryProblem(const std::string& configPath);
  virtual ~TrajectoryProblem();

  std::string getCstrName(const size_t i) const;

  // getters
  const size_t& nBoxes() const { return nBoxes_; }
  const size_t& nPlans() const { return nPlans_; }
  const size_t& nObstacles() const { return nObstacles_; }
  const size_t& nFixedPlanes() const { return nFixedPlanes_; }
  const size_t& nFixedPlanCstr() const { return nFixedPlanCstr_; }
  const size_t& nMobilePlanCstr() const { return nMobilePlanCstr_; }
  const size_t& numberOfCstr() const { return numberOfCstr_; }
  const Box& initBox() const { return initBox_; }
  const Box& finalBox() const { return finalBox_; }
  const std::vector<Box>& boxes() const { return boxes_; }
  const std::vector<Box>& obstacles() const { return obstacles_; }
  const std::vector<FixedPlan>& fixedPlanes() const { return fixedPlanes_; }
  const std::vector<PlanForHull>& plans() const { return plans_; }
  const Eigen::Vector3d& initPos() const { return initPos_; }
  const Eigen::Vector3d& finalPos() const { return finalPos_; }
  const Eigen::Vector3d& boxSize() const { return boxSize_; }
  const std::vector<std::string>& cstrNames() const { return cstrNames_; };
  const ProblemConfig& config() const { return config_; }
  const Index& dimBox() const { return dimBox_; }
  const Index& dimPlan() const { return dimPlan_; }
  const Index& dimNormal() const { return dimNormal_; }
  const Index& dimDistance() const { return dimDistance_; }
  const Index& dimBoxes() const { return dimBoxes_; }
  const Index& dimPlans() const { return dimPlans_; }
  const Index& dimNormals() const { return dimNormals_; }
  const Index& dimDistances() const { return dimDistances_; }
  const Index& dimVar() const { return dimVar_; }
  const size_t& maxIter() const { return maxIter_; }
  const double& securityDistance() const { return securityDistance_; }
  const CostDistance& costFct() const { return costFct_; };
  Eigen::VectorXd findInitPoint() const;

  const Box& getBox(Index i) const
  {
    if (i == -1)
      return initBox_;
    else if (i == nBoxes_)
      return finalBox_;
    else if (i > -1 && i < nBoxes_)
      return boxes().at(i);
    else
      std::cerr << "No box with index " << i << std::endl;
  }

  const std::vector<BoxAboveFixedPlan>& boxAboveFixedPlanFcts() const
  {
    return boxAboveFixedPlanFcts_;
  }
  const BoxAbovePlan& initBoxAbovePlanFct() const
  {
    return initBoxAbovePlanFct_;
  }
  const BoxAbovePlan& finalBoxAbovePlanFct() const
  {
    return finalBoxAbovePlanFct_;
  }
  const std::vector<BoxAbovePlan>& boxAbovePlanFcts() const
  {
    return boxAbovePlanFcts_;
  }
  const std::vector<BoxAbovePlan>& obstacleAbovePlanFcts() const
  {
    return obstacleAbovePlanFcts_;
  }

  Eigen::VectorXd getBoxPositionsFromX(const Eigen::VectorXd& x) const;
  Eigen::VectorXd getPlansNormalsFromX(const Eigen::VectorXd& x) const;
  Eigen::VectorXd getPlansDistancesFromX(const Eigen::VectorXd& x) const;
  Eigen::Vector3d getBoxPositionFromX(size_t i, const Eigen::VectorXd& x) const;
  Eigen::Vector3d getPlanNormalFromX(size_t i, const Eigen::VectorXd& x) const;
  double getPlanDistanceFromX(size_t i, const Eigen::VectorXd& x) const;

  void normalizeNormals(Eigen::Ref<Eigen::VectorXd> x) const;

  /// \brief Print method.
  /// \param o output stream.
  /// \return output stream.
  std::ostream& print (std::ostream& o) const;

 private:
  ProblemConfig config_;

  Eigen::Vector3d initPos_;
  Eigen::Vector3d finalPos_;

  std::vector<Box> obstacles_;
  std::vector<FixedPlan> fixedPlanes_;

  Eigen::Vector3d boxSize_;

  size_t nBoxes_;
  size_t nPlans_;
  size_t nObstacles_;
  size_t nFixedPlanes_;
  size_t nMobilePlanCstr_;
  size_t nFixedPlanCstr_;
  size_t numberOfCstr_;

  double securityDistance_;

  size_t maxIter_;

  Index dimBox_, dimPlan_, dimNormal_, dimDistance_;
  Index dimBoxes_, dimPlans_, dimNormals_, dimDistances_, dimVar_;

  Box initBox_, finalBox_;
  BoxAbovePlan initBoxAbovePlanFct_, finalBoxAbovePlanFct_;

  std::vector<Box> boxes_;
  std::vector<PlanForHull> plans_;
  CostDistance costFct_;
  std::vector<BoxAboveFixedPlan> boxAboveFixedPlanFcts_;
  std::vector<BoxAbovePlan> boxAbovePlanFcts_;
  std::vector<BoxAbovePlan> obstacleAbovePlanFcts_;
  std::vector<std::string> cstrNames_;
};

/// \brief Output stream operator for frames.
/// \param o output stream.
/// \param f TrajectoryProblem.
/// \return output stream.
std::ostream& operator<<(std::ostream& o, const TrajectoryProblem& f);

} /* feettrajectory */

