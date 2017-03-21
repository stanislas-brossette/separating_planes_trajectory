#include <feet-trajectory/utils/AlternateQPSolver.hh>

namespace feettrajectory
{
AlternateQPSolver::AlternateQPSolver(const TrajectoryProblem& pb)
    : pb_(pb), qpNfixed_(pb)
{
}

AlternateQPSolver::~AlternateQPSolver() {}

void AlternateQPSolver::init(const Eigen::VectorXd& xInit)
{
  qpNfixed_.formQP(pb_, pb_.getPlansNormalsFromX(xInit));
}

void AlternateQPSolver::solve() {}

} /* feettrajectory */
