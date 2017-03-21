#include <feet-trajectory/utils/QPWithNFixed.hh>

namespace feettrajectory
{
QPWithNFixed::QPWithNFixed(const TrajectoryProblem& pb)
{
  dimBox_ = 3;
  dimBoxes_ = pb.nBoxes() * dimBox_;
  dimVar_ = pb.nBoxes() * dimBox_ + pb.nPlans() * 1;
  dimCstr_ = pb.numberOfCstr();
  setDimensions(dimVar_, dimCstr_);
}
QPWithNFixed::~QPWithNFixed() {}

void QPWithNFixed::formQP(const TrajectoryProblem& pb,
                          const Eigen::VectorXd& normals)
{
  pb.costFct().fillQuadCost(A_.topLeftCorner(dimBoxes_, dimBoxes_),
                            c_.head(dimBoxes_));

  long cstrIndexBegin = 0;
  for (size_t i = 0; i < pb.boxAboveFixedPlanFcts().size(); ++i)
  {
    long boxIndexBegin = 3 * pb.boxAboveFixedPlanFcts().at(i).box().index();
    pb.boxAboveFixedPlanFcts().at(i).fillLinCstr(
        l_(cstrIndexBegin), C_.block(cstrIndexBegin, boxIndexBegin, 1, 3));
    cstrIndexBegin += 1;
  }

  for (size_t iPlan = 0; iPlan < pb.nMobilePlanCstr(); ++iPlan)
  {
    const long iBox0Above(pb.plans().at(iPlan).box0Above());
    const long iBox1Above(pb.plans().at(iPlan).box1Above());
    const long iBoxBelow(pb.plans().at(iPlan).boxBelow());
    const size_t iBox0AboveSize_t(static_cast<size_t>(iBox0Above));
    const size_t iBox1AboveSize_t(static_cast<size_t>(iBox1Above));
    const size_t iBoxBelowSize_t(static_cast<size_t>(iBoxBelow));

    Eigen::Vector3d normal;
    normal << normals.segment(3 * iPlan, 3);

    Eigen::Matrix<double, 1, 3> placeholderMatrix;

    long dIndex = dimBoxes_ + iPlan;
    long box0IndexBegin = 3 * iBox0Above;
    long box1IndexBegin = 3 * iBox1Above;

    if (iBox0Above != -1 && iBox1AboveSize_t != pb.nBoxes())
    {
      pb.boxAbovePlanFcts()
          .at(iBox0AboveSize_t)
          .fillLinCstrNFixed(l_(cstrIndexBegin),
                             C_.block(cstrIndexBegin, box0IndexBegin, 1, 3),
                             C_(cstrIndexBegin, dIndex), normal);
      cstrIndexBegin += 1;
      pb.boxAbovePlanFcts()
          .at(iBox1AboveSize_t)
          .fillLinCstrNFixed(l_(cstrIndexBegin),
                             C_.block(cstrIndexBegin, box1IndexBegin, 1, 3),
                             C_(cstrIndexBegin, dIndex), normal);
      cstrIndexBegin += 1;
    }
    else if (iBox0Above == -1 && iBox1AboveSize_t != pb.nBoxes())
    {
      pb.initBoxAbovePlanFct().fillLinCstrNFixed(
          l_(cstrIndexBegin), placeholderMatrix, C_(cstrIndexBegin, dIndex),
          normal);
      cstrIndexBegin += 1;
      pb.boxAbovePlanFcts()
          .at(iBox1AboveSize_t)
          .fillLinCstrNFixed(l_(cstrIndexBegin),
                             C_.block(cstrIndexBegin, box1IndexBegin, 1, 3),
                             C_(cstrIndexBegin, dIndex), normal);
      cstrIndexBegin += 1;
    }
    else if (iBox0Above != -1 && iBox1AboveSize_t == pb.nBoxes())
    {
      pb.boxAbovePlanFcts()
          .at(iBox0AboveSize_t)
          .fillLinCstrNFixed(l_(cstrIndexBegin),
                             C_.block(cstrIndexBegin, box0IndexBegin, 1, 3),
                             C_(cstrIndexBegin, dIndex), normal);
      cstrIndexBegin += 1;
      pb.finalBoxAbovePlanFct().fillLinCstrNFixed(
          l_(cstrIndexBegin), placeholderMatrix, C_(cstrIndexBegin, dIndex),
          normal);
      cstrIndexBegin += 1;
    }
    else
    {
      std::cerr << "Box0 is initial pos and Box1 is final, there is no "
                   "intermediary boxes..." << std::endl;
    }

    pb.obstacleAbovePlanFcts()
        .at(iBoxBelowSize_t)
        .fillLinCstrNFixed(l_(cstrIndexBegin),
                           C_.block(cstrIndexBegin, box1IndexBegin, 1, 3),
                           C_(cstrIndexBegin, dIndex), normal, true);

    cstrIndexBegin += 1;
  }
}

} /* feettrajectory */
