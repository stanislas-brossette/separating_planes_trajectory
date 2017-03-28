#include <feet-trajectory/utils/QPBoxesFixed.hh>

namespace feettrajectory
{
QPBoxesFixed::QPBoxesFixed(const TrajectoryProblem& pb) : pb_(pb)
{
  // The problem's dimension is all the variables of the planes plus the
  // relaxation term.
  // Each box above plane constraint results in 8 dim-1 constraints, thus 24 per
  // plane. And each norm-1 constraint results in one dim-1 constraint
  setDimensions(pb_.dimPlans() + 1, 24 * pb_.nPlans() + pb_.nPlans());
}
QPBoxesFixed::~QPBoxesFixed() {}

void QPBoxesFixed::addRelaxationTerm(const double& alpha)
{
  c_.tail(1) << alpha * 1;
  C_.col(dimVar_ - 1).setConstant(1);
  lVar_.tail(1) << 0;
}

Eigen::Vector3d QPBoxesFixed::getBoxPos(Index iBox, ConstRefVec xBoxes)
{
  Eigen::Vector3d boxPos;
  if (iBox == -1)
    boxPos = pb_.initPos();
  else if (iBox == pb_.nBoxes())
    boxPos = pb_.finalPos();
  else if (iBox > -1 && iBox < pb_.nBoxes())
    boxPos = xBoxes.segment(3 * iBox, 3);
  else
    std::cerr << "No box with index " << iBox << std::endl;
  return boxPos;
}

void QPBoxesFixed::formQP(ConstRefVec xBoxes, ConstRefVec xPreviousPlanes)
{
  // std::cout << "\n\nFormLP:" << std::endl;
  addRelaxationTerm(10);
  long cstrIndexBegin = 0;

  l_.setZero();

  for (size_t iPlan = 0; iPlan < pb_.nMobilePlanCstr(); ++iPlan)
  {
    Index distanceIndexBegin = 4 * iPlan;
    Index normalIndexBegin = 4 * iPlan + 1;

    const long iBox0Above(pb_.plans().at(iPlan).box0Above());
    const Box& box0 = pb_.getBox(iBox0Above);
    // std::cout << "box0: " << box0 << std::endl;
    Eigen::Vector3d box0Pos = getBoxPos(iBox0Above, xBoxes);
    // std::cout << "box0Pos: " << box0Pos << std::endl;
    for (size_t i = 0; i < box0.vertex().size(); i++)
    {
      C_.block(cstrIndexBegin, normalIndexBegin, 1, 3)
          << (box0Pos + box0.vertex().at(i)).transpose();
      C_.block(cstrIndexBegin, distanceIndexBegin, 1, 1) << -1;
      cstrIndexBegin++;
    }
    // std::cout << "C_.block: "
    //<< C_.block(cstrIndexBegin-8, normalIndexBegin, 8, 3)
    //<< std::endl;

    const long iBox1Above(pb_.plans().at(iPlan).box1Above());
    const Box& box1 = pb_.getBox(iBox1Above);
    // std::cout << "box1: " << box1 << std::endl;
    Eigen::Vector3d box1Pos = getBoxPos(iBox1Above, xBoxes);
    // std::cout << "box1Pos: " << box1Pos << std::endl;
    for (size_t i = 0; i < box1.vertex().size(); i++)
    {
      C_.block(cstrIndexBegin, normalIndexBegin, 1, 3)
          << (box1Pos + box1.vertex().at(i)).transpose();
      C_.block(cstrIndexBegin, distanceIndexBegin, 1, 1) << -1;
      cstrIndexBegin++;
    }
    // std::cout << "C_.block: "
    //<< C_.block(cstrIndexBegin-8, normalIndexBegin, 8, 3)
    //<< std::endl;

    const long iBoxBelow(pb_.plans().at(iPlan).boxBelow());
    const Box& obs = pb_.obstacles().at(iBoxBelow);
    // std::cout << "obs: " << obs << std::endl;
    Eigen::Vector3d obsPos = obs.center();
    // std::cout << "obsPos: " << obsPos << std::endl;
    for (size_t i = 0; i < obs.vertex().size(); i++)
    {
      C_.block(cstrIndexBegin, normalIndexBegin, 1, 3)
          << (-obsPos - obs.vertex().at(i)).transpose();
      C_.block(cstrIndexBegin, distanceIndexBegin, 1, 1) << 1;
      cstrIndexBegin++;
    }
    // std::cout << "C_.block: "
    //<< C_.block(cstrIndexBegin-8, normalIndexBegin, 8, 3)
    //<< std::endl;

    // Adding norm 1 constraint approximation
    C_.block(cstrIndexBegin, normalIndexBegin, 1, 3)
        << xPreviousPlanes.segment(normalIndexBegin, 3).transpose();
    C_(cstrIndexBegin, C_.cols() - 1) = 0;
    l_(cstrIndexBegin) = 0.5;
    u_(cstrIndexBegin) = 1.0;
    cstrIndexBegin++;
  }
}

} /* feettrajectory */
