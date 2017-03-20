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
  const std::vector<Box>& boxes() const { return boxes_; }
  const std::vector<Box>& obstacles() const { return obstacles_; }
  const std::vector<FixedPlan>& fixedPlanes() const { return fixedPlanes_; }
  const std::vector<PlanForHull>& plans() const { return plans_; }
  const Eigen::Vector3d& initPos() const { return initPos_; }
  const Eigen::Vector3d& finalPos() const { return finalPos_; }
  const Eigen::Vector3d& boxSize() const { return boxSize_; }
  const std::vector<std::string> cstrNames() const { return cstrNames_; };

  const std::vector<BoxAboveFixedPlan> boxAboveFixedPlanFcts() const
  {
    return boxAboveFixedPlanFcts_;
  }
  const std::vector<BoxAbovePlan> boxAbovePlanFcts() const
  {
    return boxAbovePlanFcts_;
  }
  const std::vector<BoxAbovePlan> obstacleAbovePlanFcts() const
  {
    return obstacleAbovePlanFcts_;
  }
  const CostDistance& costFct() const { return costFct_; };

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

  Box initBox_, finalBox_;

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

