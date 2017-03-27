#include <feet-trajectory/utils/AlternateQPSolver.hh>

namespace feettrajectory
{
AlternateQPSolver::AlternateQPSolver(const TrajectoryProblem& pb)
    : pb_(pb), qpPlanesFixed_(pb_), qpBoxesFixed_(pb_), res_(pb_.dimVar())
{
}

AlternateQPSolver::~AlternateQPSolver() {}

void AlternateQPSolver::init(const Eigen::VectorXd& xInit)
{
  res_ << xInit;
  qpPlanesFixed_.formQP(xInit.tail(pb_.dimPlans()));
  std::cout << "qpPlanesFixed_: " << qpPlanesFixed_ << std::endl;
  qpBoxesFixed_.formQP(xInit.head(pb_.dimBoxes()), xInit.tail(pb_.dimPlans()));
  std::cout << "qpBoxesFixed_: " << qpBoxesFixed_ << std::endl;
  QPSolver_.resize(qpPlanesFixed_.dimVar(), qpPlanesFixed_.dimCstr(),
                   Eigen::lssol::eType::QP2);
  LPSolver_.resize(qpBoxesFixed_.dimVar(), qpBoxesFixed_.dimCstr());
}

void AlternateQPSolver::formAndSolveQPPlanesFixed(RefVec x)
{
  qpPlanesFixed_.formQP(x.tail(pb_.dimPlans()));
  Eigen::MatrixXd Acopy(qpPlanesFixed_.A());
  QPSolver_.solve(qpPlanesFixed_.lVar(), qpPlanesFixed_.uVar(), Acopy,
                  qpPlanesFixed_.c(), qpPlanesFixed_.C(), qpPlanesFixed_.l(),
                  qpPlanesFixed_.u());
  if (!(QPSolver_.inform() == 0 || QPSolver_.inform() == 1))
  {
    QPSolver_.print_inform();
    std::cout << "QP solver FAILED!!! Damnit" << std::endl;
  }
  std::cout << "QPSolver_.result(): " << QPSolver_.result().transpose()
            << std::endl;
  x.head(pb_.dimBoxes()) << QPSolver_.result().head(pb_.dimBoxes());
}

void AlternateQPSolver::formAndSolveLPBoxesFixed(RefVec x)
{
  qpBoxesFixed_.formQP(x.head(pb_.dimBoxes()), x.tail(pb_.dimPlans()));
  LPSolver_.solve(qpBoxesFixed_.lVar(), qpBoxesFixed_.uVar(), qpBoxesFixed_.c(),
                  qpBoxesFixed_.C(), qpBoxesFixed_.l(), qpBoxesFixed_.u());
  if (!(LPSolver_.inform() == 0 || LPSolver_.inform() == 1))
  {
    LPSolver_.print_inform();
    std::cout << "LP solver FAILED!!! Damnit" << std::endl;
  }

  std::cout << "LPSolver_.result(): " << LPSolver_.result().transpose()
            << std::endl;
  x.tail(pb_.dimPlans()) << LPSolver_.result().head(pb_.dimPlans());
}

void AlternateQPSolver::solve()
{
  Eigen::VectorXd resFixedPlanes(qpPlanesFixed_.dimVar());
  Eigen::VectorXd resFixedBoxes(qpBoxesFixed_.dimVar());

  int nIter = 0;
  while (nIter < 10)
  {
    formAndSolveQPPlanesFixed(res_);
    std::cout << "res_.transpose(): " << res_.transpose() << std::endl;
    formAndSolveLPBoxesFixed(res_);
    std::cout << "res_.transpose(): " << res_.transpose() << std::endl;
    pb_.normalizeNormals(res_);
    nIter++;
  }
}

} /* feettrajectory */
