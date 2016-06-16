#include <iostream>

#include <Eigen/Core>

#include <cube-stacks/functions/CubeAboveFixedPlan.hh>
#include <cube-stacks/utils/quat2mat.hh>

namespace cubestacks
{
CubeAboveFixedPlan::CubeAboveFixedPlan(const Cube& c,
                                       const Eigen::Vector3d& normal,
                                       const double& d)
    : CubeAbovePlan(c), normal_(normal), d_(d)
{
}

CubeAboveFixedPlan::~CubeAboveFixedPlan() {}

void CubeAboveFixedPlan::compute(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res,
                                 const Eigen::Ref<const Eigen::Vector3d> trans,
                                 const Eigen::Ref<const Eigen::Vector4d> quat) const
{
  CubeAbovePlan::compute(res, trans, quat, d_, normal_);
}

void CubeAboveFixedPlan::diffTrans(
    Eigen::Ref<Eigen::Matrix<double, 8, 3>> res,
    const Eigen::Ref<const Eigen::Vector3d> trans,
    const Eigen::Ref<const Eigen::Vector4d> quat) const
{
  CubeAbovePlan::diffTrans(res, trans, quat, d_, normal_);
}

void CubeAboveFixedPlan::diffQuat(Eigen::Ref<Eigen::Matrix<double, 8, 4>> res,
                                  const Eigen::Ref<const Eigen::Vector3d> trans,
                                  const Eigen::Ref<const Eigen::Vector4d> quat) const
{
  CubeAbovePlan::diffQuat(res, trans, quat, d_, normal_);
}

} /* cubestacks */
