#pragma once
#include <iostream>

#include <feet-trajectory/utils/QP.hh>
#include <feet-trajectory/utils/defs.hh>
#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/functions/PlanBetweenBoxAndObstacle.hh>

namespace feettrajectory
{
class QPBoxesFixedIndividual : public QP
{
 public:
  QPBoxesFixedIndividual(const TrajectoryProblem& pb);
  virtual ~QPBoxesFixedIndividual();
  void addRelaxationTerm(const double& alpha);
  void formQP(const size_t& iPlan, ConstRefVec boxesVariables, ConstRefVec xPreviousPlanes);

 private:
  const TrajectoryProblem& pb_;
  Eigen::Vector3d getBoxPos(Index iBox, ConstRefVec xBoxes);
};

} /* feettrajectory */

