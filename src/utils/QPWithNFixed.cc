#include <feet-trajectory/utils/QPWithNFixed.hh>

namespace feettrajectory
{
QPWithNFixed::QPWithNFixed(const TrajectoryProblem& pb)
{
  size_t dimBox = 3;
  size_t dimBoxes = pb.nBoxes() * dimBox;
  size_t dimVar = pb.nBoxes() * dimBox + pb.nPlans() * 1;
  size_t dimCstr = pb.boxAbovePlanFcts().size() * 8 +
                   pb.boxAboveFixedPlanFcts().size() * 8 +
                   pb.obstacleAbovePlanFcts().size() * 8;
  setDimensions(dimVar, dimCstr);
  A_.block(0, 0, dimBoxes, dimBoxes) = pb.costFct().Q();
  c_.head(dimBoxes) = pb.costFct().c();
  std::cout << "QPwithNFixed_.A_: " << A_ << std::endl;
  std::cout << "QPwithNFixed_.c_: " << c_ << std::endl;

  long cstrIndexBegin = 0;
  for (size_t i = 0; i < pb.boxAboveFixedPlanFcts().size(); ++i)
  {
    long boxIndexBegin = 3 * pb.boxAboveFixedPlanFcts().at(i).box().index();
    pb.boxAboveFixedPlanFcts().at(i).fillLinCstr(
        l_.segment(cstrIndexBegin, 8),
        C_.block(cstrIndexBegin, boxIndexBegin, 8, 3));
    cstrIndexBegin += 8;
  }
  std::cout << "l_: " << l_ << std::endl;
  std::cout << "C_: " << C_ << std::endl;
  // for (size_t i = 0; i < pb.boxAbovePlanFcts(); ++i)
  //{
  ////if (pb.boxAbovePlanFcts().box().fixed())
  ////pb.boxAbovePlanFcts().fillLinCstrNFixedBoxFixed(
  ////QPwithNFixed_.l_.segment(cstrIndexBegin, 8),
  ////C_.block(3 * pb.nBoxes() + ));
  ////cstrIndexBegin += 8;
  //}
}
QPWithNFixed::~QPWithNFixed() {}
} /* feettrajectory */
