#include <feet-trajectory/utils/AlternateQPSolverJerk.hh>
#include <feet-trajectory/utils/defs.hh>

namespace feettrajectory
{
AlternateQPSolverJerk::AlternateQPSolverJerk(const TrajectoryProblem& pb,
                                             const size_t& maxIter,
                                             const Jerk3DIntegrator& integ,
                                             const Eigen::VectorXd& state0)
    : pb_(pb),
      integ_(integ),
      state0_(state0),
      qpPlanesFixed_(pb_),
      qpBoxesFixed_(pb_),
      qpBoxesFixedIndividual_(pb_),
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
  qpPlanesFixed_.formQP(x.tail(pb_.dimPlans()));
  qpPlanesFixed_.updatePlanD(x.tail(pb_.dimPlans()));
  // TODO: Modify qp on pos to qp on jerk

  Eigen::MatrixXd Acopy(qpPlanesFixed_.A());
  Acopy.block(0, 0, pb_.dimBoxes(), pb_.dimBoxes()).setIdentity();
  Eigen::MatrixXd CJerk(qpPlanesFixed_.C());
  Eigen::VectorXd lJerk(qpPlanesFixed_.l());
  Eigen::VectorXd uJerk(qpPlanesFixed_.u());

  CJerk.leftCols(pb_.dimBoxes()) = qpPlanesFixed_.C().leftCols(pb_.dimBoxes()) *
                                   integ_.SelPos() * integ_.Uu();

  lJerk.head(qpPlanesFixed_.dimCstr()) =
      qpPlanesFixed_.l().head(qpPlanesFixed_.dimCstr()) -
      qpPlanesFixed_.C().leftCols(pb_.dimBoxes()) * integ_.SelPos() *
          integ_.Ux() * state0_;

  uJerk.head(qpPlanesFixed_.dimCstr()) =
      qpPlanesFixed_.u().head(qpPlanesFixed_.dimCstr()) -
      qpPlanesFixed_.C().leftCols(pb_.dimBoxes()) * integ_.SelPos() *
          integ_.Ux() * state0_;

  QPSolver_.solve(qpPlanesFixed_.lVar(), qpPlanesFixed_.uVar(), Acopy,
                  qpPlanesFixed_.c(), CJerk, lJerk, uJerk);
  if (!(QPSolver_.inform() == 0 || QPSolver_.inform() == 1))
  {
    QPSolver_.print_inform();
    std::cerr << "QP solver FAILED!!! Damnit" << std::endl;
  }
  // std::cout << "QPSolver_.result(): \n" << QPSolver_.result().transpose()
  //<< std::endl;
  x.head(pb_.dimBoxes()) << QPSolver_.result().head(pb_.dimBoxes());
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
    // std::cout << "LPSolverIndiv_.result(): \n"
    //<< LPSolverIndiv_.result().transpose().format(fmt::custom)
    //<< std::endl;
    x.segment(pb_.dimBoxes() + 4 * iPlan, 4) << LPSolverIndiv_.result().head(4);
  }
}

void AlternateQPSolverJerk::solve()
{
  Eigen::VectorXd resFixedPlanes(qpPlanesFixed_.dimVar());
  Eigen::VectorXd resFixedBoxes(qpBoxesFixed_.dimVar());

  int nIter = 1;
  while (nIter < maxIter_)
  {
    formAndSolveIndividualLPBoxesFixed(res_);
    //{
    // Eigen::VectorXd copyRes(res_);
    // formAndSolveLPBoxesFixed(copyRes);
    //}
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
