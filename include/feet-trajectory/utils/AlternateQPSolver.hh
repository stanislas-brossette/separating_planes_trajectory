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
  void init(const Eigen::VectorXd& xInit);
  void solve();
  const QPWithNFixed& qpNfixed() const { return qpNfixed_; }

 private:
  const TrajectoryProblem& pb_;
  QPWithNFixed qpNfixed_;
  QP qpBfixed_;

};
} /* feettrajectory */
