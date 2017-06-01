#include <feet-trajectory/functions/FixedBoxPosition.hh>
namespace feettrajectory
{
FixedBoxPosition::FixedBoxPosition(const Eigen::Vector3d& targetPos)
    : targetPos_(targetPos)
{
}
FixedBoxPosition::~FixedBoxPosition() {}

void FixedBoxPosition::fillLinCstr(const Eigen::Vector3d& targetPos, RefVec lb, RefMat C, RefVec ub)
{
  C.setIdentity();
  lb << targetPos;
  ub << targetPos;
}

void FixedBoxPosition::compute(
    Eigen::Ref<Eigen::Matrix<double, 3, 1>> res,
    const Eigen::Ref<const Eigen::Vector3d> pos) const
{
  res = pos;
}
void FixedBoxPosition::diff(Eigen::Ref<Eigen::Matrix<double, 3, 3>> res) const
{
  res = Eigen::Matrix3d::Identity();
}
void FixedBoxPosition::LB(Eigen::Ref<Eigen::Matrix<double, 3, 1>> res) const
{
  res = targetPos_;
}
void FixedBoxPosition::UB(Eigen::Ref<Eigen::Matrix<double, 3, 1>> res) const
{
  res = targetPos_;
}
} /* feettrajectory */
