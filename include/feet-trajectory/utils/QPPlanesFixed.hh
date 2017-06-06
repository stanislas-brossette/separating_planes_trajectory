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
  //void addRelaxationTerm(const double& alpha);
  void addRelaxationTerm(size_t cstrIndexBegin, size_t cstrSize, bool isVirtual);
  void formQP(ConstRefVec planVariables);
  void updatePlanD(RefVec planVariables);

 private:
  const TrajectoryProblem& pb_;
  /// @brief relaxation term for non virtual collision avoidance
  double alpha_ = 1000;
  /// @brief relaxation term for virtual collision avoidance
  double alphaVirtual_ = 100;
};

} /* feettrajectory */

