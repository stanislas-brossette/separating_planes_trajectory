#pragma once
#include <iostream>
#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/utils/QP.hh>
#include <feet-trajectory/utils/QPWithNFixed.hh>

namespace feettrajectory
{
class AlternateQPSolver
{
 public:
  AlternateQPSolver(const TrajectoryProblem& pb);
  virtual ~AlternateQPSolver();
  void init();
  void solve();

 private:
  const TrajectoryProblem& pb_;
  QPWithNFixed QPwithNFixed_;
  QP QPwithBfixed_;

};
} /* feettrajectory */
