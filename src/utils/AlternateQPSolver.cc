#include <feet-trajectory/utils/AlternateQPSolver.hh>
#include <feet-trajectory/utils/defs.hh>

namespace feettrajectory
{
AlternateQPSolver::AlternateQPSolver(const TrajectoryProblem& pb,
                                     const size_t& maxIter)
    : pb_(pb),
      qpPlanesFixed_(pb_),
      qpBoxesFixed_(pb_),
      res_(pb_.dimVar()),
      maxIter_(maxIter)
{
  resHistory_.resize(maxIter_ + 1);
  for (size_t i = 0; i < resHistory_.size(); i++)
  {
    resHistory_[i].resize(pb_.dimVar());
    resHistory_[i].setZero();
  }
}

AlternateQPSolver::~AlternateQPSolver() {}

void AlternateQPSolver::init(const Eigen::VectorXd& xInit)
{
  res_ << xInit;
  resHistory_[0] << res_;
  qpPlanesFixed_.formQP(xInit.tail(pb_.dimPlans()));
  qpBoxesFixed_.formQP(xInit.head(pb_.dimBoxes()), xInit.tail(pb_.dimPlans()));
  QPSolver_.resize(qpPlanesFixed_.dimVar(), qpPlanesFixed_.dimCstr(),
                   Eigen::lssol::eType::QP2);
  LPSolver_.resize(qpBoxesFixed_.dimVar(), qpBoxesFixed_.dimCstr());
}

void AlternateQPSolver::formAndSolveQPPlanesFixed(RefVec x)
{
  qpPlanesFixed_.formQP(x.tail(pb_.dimPlans()));
  qpPlanesFixed_.updatePlanD(x.tail(pb_.dimPlans()));
  Eigen::MatrixXd Acopy(qpPlanesFixed_.A());
  QPSolver_.solve(qpPlanesFixed_.lVar(), qpPlanesFixed_.uVar(), Acopy,
                  qpPlanesFixed_.c(), qpPlanesFixed_.C(), qpPlanesFixed_.l(),
                  qpPlanesFixed_.u());
  if (!(QPSolver_.inform() == 0 || QPSolver_.inform() == 1))
  {
    QPSolver_.print_inform();
    std::cerr << "QP solver FAILED!!! Damnit" << std::endl;
  }
  // std::cout << "QPSolver_.result(): \n" << QPSolver_.result().transpose()
  //<< std::endl;
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
    std::cerr << "LP solver FAILED!!! Damnit" << std::endl;
  }
  // std::cout << "LPSolver_.result(): \n" <<
  // LPSolver_.result().transpose().format(fmt::custom)
  //<< std::endl;
  x.tail(pb_.dimPlans()) << LPSolver_.result().head(pb_.dimPlans());
}

void AlternateQPSolver::solve()
{
  Eigen::VectorXd resFixedPlanes(qpPlanesFixed_.dimVar());
  Eigen::VectorXd resFixedBoxes(qpBoxesFixed_.dimVar());

  int nIter = 1;
  while (nIter < maxIter_)
  {
    formAndSolveLPBoxesFixed(res_);
    // std::cout << "res_: \n" << res_.transpose().format(fmt::custom) <<
    // std::endl;
    pb_.normalizeNormals(res_);
    // std::cout << "normalized res_: \n" <<
    // res_.transpose().format(fmt::custom) << std::endl;
    resHistory_[nIter] << res_;
    nIter++;

    // std::cout << "Solving QP planes fixed" << std::endl;
    formAndSolveQPPlanesFixed(res_);
    // std::cout << "res_: \n" << res_.transpose() << std::endl;
    resHistory_[nIter] << res_;
    nIter++;
  }
  totalIter_ = nIter;
}

void AlternateQPSolver::logAllX(const std::string& folderName) const
{
  std::ofstream xLogFile;
  xLogFile.open(folderName + "xLog.m");

  for (size_t i = 0; i < totalIter_; i++)
  {
    xLogFile << "%============== iteration " << i
             << "==================" << std::endl;
    xLogFile << "x_" << i << " = ";
    xLogFile << resHistory_[i].format(fmt::matlab) << std::endl;
  }
  xLogFile.close();
}

} /* feettrajectory */
