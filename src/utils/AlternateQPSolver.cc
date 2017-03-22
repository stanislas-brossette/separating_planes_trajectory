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
  QPSolver_.resize(pb_.dimVar(), pb_.numberOfCstr(), Eigen::lssol::eType::QP2);
}

void AlternateQPSolver::solve()
{
  Eigen::VectorXd res(qpNfixed_.dimVar());
  Eigen::MatrixXd Acopy(qpNfixed_.A());
  QPSolver_.print(qpNfixed_.lVar(), qpNfixed_.uVar(), Acopy, qpNfixed_.c(),
                  qpNfixed_.C(), qpNfixed_.l(), qpNfixed_.u());
  QPSolver_.solve(qpNfixed_.lVar(), qpNfixed_.uVar(), Acopy, qpNfixed_.c(),
                  qpNfixed_.C(), qpNfixed_.l(), qpNfixed_.u());
  if (!(QPSolver_.inform() == 0 || QPSolver_.inform() == 1))
  {
    QPSolver_.print_inform();
    std::cout << "QP solver FAILED!!! Damnit" << std::endl;
  }
  else
  {
    res = QPSolver_.result();
    std::cout << "res: " << res << std::endl;
  }
}

} /* feettrajectory */
