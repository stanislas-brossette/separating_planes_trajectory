#include <feet-trajectory/functions/FixedBoxPosition.hh>
namespace feettrajectory
{
FixedBoxPosition::FixedBoxPosition(const Box& box,
                                   const Eigen::Vector3d& targetPos)
    : box_(box), targetPos_(targetPos)
{
}
FixedBoxPosition::~FixedBoxPosition() {}

void FixedBoxPosition::fillLinCstr(const Eigen::Vector3d& targetPos, RefVec lb, RefMat C, RefVec ub)
{
  C.setIdentity();
  lb << targetPos;
  ub << targetPos;
}
} /* feettrajectory */
