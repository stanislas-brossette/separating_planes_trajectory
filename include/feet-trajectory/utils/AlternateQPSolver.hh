#pragma once

#include <iostream>
#include <fstream>

#include <eigenqp/LSSOL.h>
#include <eigenqp/LSSOL_QP.h>
#include <eigenqp/LSSOL_LP.h>

#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/utils/defs.hh>
#include <feet-trajectory/utils/QP.hh>
#include <feet-trajectory/utils/QPPlanesFixed.hh>
#include <feet-trajectory/utils/QPBoxesFixed.hh>
#include <feet-trajectory/utils/QPBoxesFixedIndividual.hh>

namespace feettrajectory
{
class AlternateQPSolver
{
 public:
  AlternateQPSolver(const TrajectoryProblem& pb, const size_t& maxIter, const double& precision = 1e-8);
  virtual ~AlternateQPSolver();
  void init(const Eigen::VectorXd& xInit);
  void solve();
  void logAllX(const std::string& fileName) const;
  const QPPlanesFixed& qpPlanesFixed() const { return qpPlanesFixed_; }
  const QPBoxesFixed& qpBoxesFixed() const { return qpBoxesFixed_; }
  void formAndSolveQPPlanesFixed(RefVec x);
  void formAndSolveLPBoxesFixed(RefVec x);
  void formAndSolveIndividualLPBoxesFixed(RefVec x);
  Eigen::VectorXd res() const { return res_; }
  const size_t& totalIter() const { return totalIter_; }

 private:
  const TrajectoryProblem& pb_;
  QPPlanesFixed qpPlanesFixed_;
  QPBoxesFixed qpBoxesFixed_;
  QPBoxesFixedIndividual qpBoxesFixedIndividual_;
  Eigen::VectorXd res_;
  size_t maxIter_, totalIter_;
  std::vector<Eigen::VectorXd> resHistory_;
  double precision_;

  /// @brief QP solver
  Eigen::LSSOL_QP QPSolver_;
  Eigen::LSSOL_LP LPSolver_;
  Eigen::LSSOL_LP LPSolverIndiv_;

};
} /* feettrajectory */
