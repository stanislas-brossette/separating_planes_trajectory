#include <feet-trajectory/utils/AlternateQPSolverJerk.hh>
#include <feet-trajectory/utils/defs.hh>

namespace feettrajectory
{
AlternateQPSolverJerk::AlternateQPSolverJerk(const TrajectoryProblem& pb,
                                             const size_t& maxIter,
                                             const Jerk3DIntegrator& integ,
                                             const Eigen::VectorXd& state0,
                                             const double& prec)
    : pb_(pb),
      integ_(integ),
      state0_(state0),
      qpPlanesFixed_(pb_),
      qpBoxesFixed_(pb_),
      qpBoxesFixedIndividual_(pb_),
      res_(pb_.dimVar()),
      maxIter_(maxIter),
      precision_(prec)
{
  resHistory_.resize(maxIter_ + 1);
  for (size_t i = 0; i < resHistory_.size(); i++)
  {
    resHistory_[i].resize(pb_.dimVar());
    resHistory_[i].setZero();
  }
}

AlternateQPSolverJerk::~AlternateQPSolverJerk() {}

void AlternateQPSolverJerk::init(const Eigen::VectorXd& xInit)
{
  res_ << xInit;
  resHistory_[0] << res_;
  qpPlanesFixed_.formQP(xInit.tail(pb_.dimPlans()));
  qpBoxesFixed_.formQP(integ_.getPos(xInit.head(pb_.dimBoxes()), state0_),
                       xInit.tail(pb_.dimPlans()));
  QPSolver_.resize(qpPlanesFixed_.dimVar(), qpPlanesFixed_.dimCstr(),
                   Eigen::lssol::eType::QP2);
  LPSolver_.resize(qpBoxesFixed_.dimVar(), qpBoxesFixed_.dimCstr());
  LPSolverIndiv_.resize(5, 25);
}

void AlternateQPSolverJerk::formAndSolveQPPlanesFixed(RefVec x)
{
  auto dimBoxes = pb_.dimBoxes();

  qpPlanesFixed_.formQP(x.tail(pb_.dimPlans()));
  qpPlanesFixed_.updatePlanD(x.tail(pb_.dimPlans()));

  Eigen::MatrixXd AJerk(qpPlanesFixed_.A());
  Eigen::MatrixXd CJerk(qpPlanesFixed_.C());
  Eigen::VectorXd lJerk(qpPlanesFixed_.l());
  Eigen::VectorXd uJerk(qpPlanesFixed_.u());
  Eigen::VectorXd cJerk(qpPlanesFixed_.c());

  Eigen::MatrixXd SposUu = integ_.SelPos()*integ_.Uu();
  Eigen::MatrixXd SposUx = integ_.SelPos()*integ_.Ux();

  AJerk.block(0, 0, dimBoxes, dimBoxes) =
      Eigen::MatrixXd::Identity(dimBoxes, dimBoxes) +
      SposUu.transpose() * qpPlanesFixed_.A().block(0, 0, dimBoxes, dimBoxes) *
          SposUu;
  cJerk.head(dimBoxes) =
      SposUu.transpose() * qpPlanesFixed_.c().head(dimBoxes) +
      SposUu.transpose() * qpPlanesFixed_.A().block(0, 0, dimBoxes, dimBoxes) *
          SposUx * state0_;

  CJerk.leftCols(dimBoxes) =
      qpPlanesFixed_.C().leftCols(dimBoxes) * SposUu;

  lJerk = qpPlanesFixed_.l() -
          qpPlanesFixed_.C().leftCols(dimBoxes) * SposUx * state0_;

  uJerk = qpPlanesFixed_.u() -
          qpPlanesFixed_.C().leftCols(dimBoxes) * SposUx * state0_;

  QPSolver_.solve(qpPlanesFixed_.lVar(), qpPlanesFixed_.uVar(), AJerk,
                  cJerk, CJerk, lJerk, uJerk);
  if (!(QPSolver_.inform() == 0 || QPSolver_.inform() == 1))
  {
    QPSolver_.print_inform();
    std::cerr << "QP solver FAILED!!! Damnit" << std::endl;
  }
  // std::cout << "QPSolver_.result(): \n" << QPSolver_.result().transpose()
  //<< std::endl;
  x.head(dimBoxes) << QPSolver_.result().head(dimBoxes);
}

void AlternateQPSolverJerk::formAndSolveLPBoxesFixed(RefVec x)
{
  qpBoxesFixed_.formQP(integ_.getPos(x.head(pb_.dimBoxes()), state0_),
                       x.tail(pb_.dimPlans()));
  LPSolver_.solve(qpBoxesFixed_.lVar(), qpBoxesFixed_.uVar(), qpBoxesFixed_.c(),
                  qpBoxesFixed_.C(), qpBoxesFixed_.l(), qpBoxesFixed_.u());
  if (!(LPSolver_.inform() == 0 || LPSolver_.inform() == 1))
  {
    LPSolver_.print_inform();
    std::cerr << "LP solver FAILED!!! Damnit" << std::endl;
  }
  // std::cout << "LPSolver_.result(): \n"
  //<< LPSolver_.result().transpose().format(fmt::custom) << std::endl;
  // x.tail(pb_.dimPlans()) << LPSolver_.result().head(pb_.dimPlans());
}

void AlternateQPSolverJerk::formAndSolveIndividualLPBoxesFixed(RefVec x)
{
  for (size_t iPlan = 0; iPlan < pb_.nMobilePlanCstr(); ++iPlan)
  {
    qpBoxesFixedIndividual_.formQP(
        iPlan, integ_.getPos(x.head(pb_.dimBoxes()), state0_),
        x.tail(pb_.dimPlans()));
    // std::cout << "qpBoxesFixedIndividual_ for plan " << iPlan << ": ";
    // qpBoxesFixedIndividual_.print(std::cout, fmt::matlab);
    LPSolverIndiv_.solve(
        qpBoxesFixedIndividual_.lVar(), qpBoxesFixedIndividual_.uVar(),
        qpBoxesFixedIndividual_.c(), qpBoxesFixedIndividual_.C(),
        qpBoxesFixedIndividual_.l(), qpBoxesFixedIndividual_.u());
    if (!(LPSolverIndiv_.inform() == 0 || LPSolverIndiv_.inform() == 1))
    {
      LPSolverIndiv_.print_inform();
      std::cerr << "LP solver FAILED!!! Damnit" << std::endl;
    }
    x.segment(pb_.dimBoxes() + 4 * iPlan, 4) << LPSolverIndiv_.result().head(4);
  }
}

void AlternateQPSolverJerk::solve()
{
  Eigen::VectorXd prevRes(res_.rows());
  prevRes.setZero();
  bool converged = false;

  int nIter = 1;
  while (nIter < maxIter_ && !converged)
  {
    formAndSolveIndividualLPBoxesFixed(res_);
    pb_.normalizeNormals(res_);
    resHistory_[nIter] << res_;
    nIter++;

    formAndSolveQPPlanesFixed(res_);
    resHistory_[nIter] << res_;

    if((prevRes - res_).head(pb_.dimBoxes()).norm() < precision_)
    {
      converged = true;
      std::cout << "Alternate QP Jerk CONVERGED on iteration " << nIter << std::endl;
    }
    prevRes = res_;

    nIter++;
  }
  totalIter_ = nIter;
}

void AlternateQPSolverJerk::logAllX(const std::string& folderName) const
{
  std::ofstream xLogFile;
  xLogFile.open(folderName + "xLog.m");

  Eigen::VectorXd posAtResi(resHistory_[0]);
  for (size_t i = 0; i < totalIter_; i++)
  {
    xLogFile << "%============== iteration " << i
             << "==================" << std::endl;
    xLogFile << "x_" << i << " = ";
    posAtResi = resHistory_[i];
    posAtResi.head(pb_.dimBoxes()) =
        integ_.getPos(resHistory_[i].head(pb_.dimBoxes()), state0_);
    xLogFile << posAtResi.format(fmt::matlab) << std::endl;
  }
  xLogFile.close();
}

} /* feettrajectory */
