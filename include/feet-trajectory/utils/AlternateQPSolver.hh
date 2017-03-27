#pragma once

#include <iostream>
#include <fstream>

#include <EigenQP/LSSOL.h>
#include <EigenQP/LSSOL_QP.h>
#include <EigenQP/LSSOL_LP.h>

#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/utils/defs.hh>
#include <feet-trajectory/utils/QP.hh>
#include <feet-trajectory/utils/QPPlanesFixed.hh>
#include <feet-trajectory/utils/QPBoxesFixed.hh>

namespace feettrajectory
{
class AlternateQPSolver
{
 public:
  AlternateQPSolver(const TrajectoryProblem& pb, const size_t& maxIter);
  virtual ~AlternateQPSolver();
  void init(const Eigen::VectorXd& xInit);
  void solve();
  void logAllX(const std::string& fileName) const;
  const QPPlanesFixed& qpPlanesFixed() const { return qpPlanesFixed_; }
  const QPBoxesFixed& qpBoxesFixed() const { return qpBoxesFixed_; }
  void formAndSolveQPPlanesFixed(RefVec x);
  void formAndSolveLPBoxesFixed(RefVec x);
  Eigen::VectorXd res() const { return res_; }

 private:
  const TrajectoryProblem& pb_;
  QPPlanesFixed qpPlanesFixed_;
  QPBoxesFixed qpBoxesFixed_;
  Eigen::VectorXd res_;
  size_t maxIter_, totalIter_;
  std::vector<Eigen::VectorXd> resHistory_;

  /// @brief QP solver
  Eigen::LSSOL_QP QPSolver_;
  Eigen::LSSOL_LP LPSolver_;

};
} /* feettrajectory */
