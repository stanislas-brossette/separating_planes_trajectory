#include <feet-trajectory/utils/QPPlanesFixed.hh>

namespace feettrajectory
{
QPPlanesFixed::QPPlanesFixed(const TrajectoryProblem& pb) : pb_(pb)
{
  setDimensions(pb_.dimBoxes() + 1,
                pb_.nFixedPlanes() * pb_.nBoxes() + 2 * pb_.nPlans());
}
QPPlanesFixed::~QPPlanesFixed() {}

void QPPlanesFixed::addRelaxationTerm(const double& alpha)
{
  c_.tail(1) << alpha * 1;
  C_.col(dimVar_ - 1).setConstant(1);
  lVar_.tail(1) << 0;
}

void QPPlanesFixed::formQP(ConstRefVec xPlanes)
{
  pb_.costFct().fillQuadCost(A_.topLeftCorner(pb_.dimBoxes(), pb_.dimBoxes()),
                             c_.head(pb_.dimBoxes()));
  A_.bottomRightCorner(1, 1) << -1;

  long cstrIndexBegin = 0;
  for (size_t i = 0; i < pb_.boxAboveFixedPlanFcts().size(); ++i)
  {
    long boxIndexBegin = 3 * pb_.boxAboveFixedPlanFcts().at(i).box().index();
    pb_.boxAboveFixedPlanFcts().at(i).fillLinCstr(
        l_(cstrIndexBegin), C_.block(cstrIndexBegin, boxIndexBegin, 1, 3));
    cstrIndexBegin += 1;
  }

  for (size_t iPlan = 0; iPlan < pb_.nMobilePlanCstr(); ++iPlan)
  {
    const long iBox0Above(pb_.plans().at(iPlan).box0Above());
    const long iBox1Above(pb_.plans().at(iPlan).box1Above());
    const long iBoxBelow(pb_.plans().at(iPlan).boxBelow());
    const size_t iBox0AboveSize_t(static_cast<size_t>(iBox0Above));
    const size_t iBox1AboveSize_t(static_cast<size_t>(iBox1Above));
    const size_t iBoxBelowSize_t(static_cast<size_t>(iBoxBelow));

    Eigen::Vector3d planN;
    double planD(xPlanes(4 * iPlan));
    planN << xPlanes.segment(4 * iPlan + 1, 3);

    Eigen::Matrix<double, 1, 3> placeholderMatrix;

    Index box0IndexBegin = 3 * iBox0Above;
    Index box1IndexBegin = 3 * iBox1Above;

    if (!pb_.getBox(iBox0Above).fixed())
    {
      PlanBetweenBoxAndObstacle::fillLinCstr(
          pb_.getBox(iBox0Above), pb_.obstacles().at(iBoxBelow), planD, planN,
          l_.segment(cstrIndexBegin, 1),
          C_.block(cstrIndexBegin, box0IndexBegin, 1, 3),
          pb_.securityDistance());
    }
    else
    {
      PlanBetweenBoxAndObstacle::fillLinCstr(
          pb_.getBox(iBox0Above), pb_.obstacles().at(iBoxBelow), planD, planN,
          l_.segment(cstrIndexBegin, 1), placeholderMatrix,
          pb_.securityDistance());
    }
    cstrIndexBegin += 1;

    if (!pb_.getBox(iBox1Above).fixed())
    {
      PlanBetweenBoxAndObstacle::fillLinCstr(
          pb_.getBox(iBox1Above), pb_.obstacles().at(iBoxBelow), planD, planN,
          l_.segment(cstrIndexBegin, 1),
          C_.block(cstrIndexBegin, box1IndexBegin, 1, 3),
          pb_.securityDistance());
    }
    else
    {
      PlanBetweenBoxAndObstacle::fillLinCstr(
          pb_.getBox(iBox1Above), pb_.obstacles().at(iBoxBelow), planD, planN,
          l_.segment(cstrIndexBegin, 1), placeholderMatrix,
          pb_.securityDistance());
    }
    cstrIndexBegin += 1;
  }

  addRelaxationTerm(10);
}

void QPPlanesFixed::updatePlanD(RefVec xPlanes)
{
  for (size_t iPlan = 0; iPlan < pb_.nMobilePlanCstr(); ++iPlan)
  {
    const long iBoxBelow(pb_.plans().at(iPlan).boxBelow());
    PlanBetweenBoxAndObstacle::updatePlanD(
        pb_.obstacles().at(iBoxBelow), xPlanes(4 * iPlan),
        xPlanes.segment(4 * iPlan + 1, 3), pb_.securityDistance());
  }
}

} /* feettrajectory */
