#pragma once

#include <iostream>
#include <feet-trajectory/utils/QP.hh>
#include <feet-trajectory/utils/defs.hh>
#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/functions/PlanBetweenBoxAndObstacle.hh>

namespace feettrajectory
{
class QPPlanesFixed : public QP
{
 public:
  QPPlanesFixed(const TrajectoryProblem& pb);
  virtual ~QPPlanesFixed();
  void addRelaxationTerm(const double& alpha);
  void formQP(ConstRefVec planVariables);

 private:
  const TrajectoryProblem& pb_;
};

} /* feettrajectory */

