#include <feet-trajectory/utils/AlternateQPSolver.hh>

namespace feettrajectory
{
AlternateQPSolver::AlternateQPSolver(const TrajectoryProblem& pb)
    : pb_(pb),
      // qpNfixed_(pb),
      qpPlanesFixed_(pb)
{
}

AlternateQPSolver::~AlternateQPSolver() {}

void AlternateQPSolver::init(const Eigen::VectorXd& xInit)
{
  // qpNfixed_.formQP(pb_, pb_.getPlansNormalsFromX(xInit));
  qpPlanesFixed_.formQP(xInit.tail(pb_.dimPlans()));
  std::cout << "qpPlanesFixed_: " << qpPlanesFixed_ << std::endl;
  QPSolver_.resize(pb_.dimVar(), pb_.numberOfCstr(), Eigen::lssol::eType::QP2);
}

void AlternateQPSolver::solve()
{
  Eigen::VectorXd resFixedPlanes(qpPlanesFixed_.dimVar());
  Eigen::MatrixXd Acopy(qpPlanesFixed_.A());
  QPSolver_.print(qpPlanesFixed_.lVar(), qpPlanesFixed_.uVar(),
                  qpPlanesFixed_.A(), qpPlanesFixed_.c(), qpPlanesFixed_.C(),
                  qpPlanesFixed_.l(), qpPlanesFixed_.u());
  QPSolver_.solve(qpPlanesFixed_.lVar(), qpPlanesFixed_.uVar(), Acopy,
                  qpPlanesFixed_.c(), qpPlanesFixed_.C(), qpPlanesFixed_.l(),
                  qpPlanesFixed_.u());
  if (!(QPSolver_.inform() == 0 || QPSolver_.inform() == 1))
  {
    QPSolver_.print_inform();
    std::cout << "QP solver FAILED!!! Damnit" << std::endl;
  }
  else
  {
    resFixedPlanes = QPSolver_.result();
    std::cout << "resFixedPlanes: " << resFixedPlanes << std::endl;
  }
}

} /* feettrajectory */
