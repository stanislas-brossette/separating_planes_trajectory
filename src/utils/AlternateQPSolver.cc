#include <feet-trajectory/utils/AlternateQPSolver.hh>

namespace feettrajectory
{
AlternateQPSolver::AlternateQPSolver(const TrajectoryProblem& pb)
    : pb_(pb), QPwithNFixed_(pb)
{
}

AlternateQPSolver::~AlternateQPSolver() {}

void AlternateQPSolver::init() {}

void AlternateQPSolver::solve() {}

} /* feettrajectory */
