#include <iostream>

#include <Eigen/Core>

#include <feet-trajectory/functions/BoxAbovePlan.hh>
#include <feet-trajectory/functions/BoxAboveFixedPlan.hh>
#include <feet-trajectory/utils/Box.hh>
#include <feet-trajectory/utils/FixedPlan.hh>
#include <feet-trajectory/utils/quat2mat.hh>

namespace feettrajectory
{
BoxAboveFixedPlan::BoxAboveFixedPlan(const Box& c,
                                     const Eigen::Vector3d& normal,
                                     const double& d)
    : BoxAbovePlan(c), normal_(normal), d_(d)
{
}
BoxAboveFixedPlan::BoxAboveFixedPlan(const Box& c, const FixedPlan& p)
    : BoxAbovePlan(c), normal_(p.normal()), d_(p.d())
{
}

BoxAboveFixedPlan::~BoxAboveFixedPlan() {}

void BoxAboveFixedPlan::compute(
    Eigen::Ref<Eigen::Matrix<double, 8, 1>> res,
    const Eigen::Ref<const Eigen::Vector3d> trans,
    const Eigen::Ref<const Eigen::Vector4d> quat) const
{
  BoxAbovePlan::compute(res, trans, quat, d_, normal_);
}

void BoxAboveFixedPlan::diffTrans(
    Eigen::Ref<Eigen::Matrix<double, 8, 3>> res,
    const Eigen::Ref<const Eigen::Vector3d> trans,
    const Eigen::Ref<const Eigen::Vector4d> quat) const
{
  BoxAbovePlan::diffTrans(res, trans, quat, d_, normal_);
}

void BoxAboveFixedPlan::diffQuat(
    Eigen::Ref<Eigen::Matrix<double, 8, 4>> res,
    const Eigen::Ref<const Eigen::Vector3d> trans,
    const Eigen::Ref<const Eigen::Vector4d> quat) const
{
  BoxAbovePlan::diffQuat(res, trans, quat, d_, normal_);
}

void BoxAboveFixedPlan::diffTrans(Eigen::Ref<Eigen::Matrix<double, 1, 3>> res,
                                  const Eigen::Ref<const Eigen::Vector3d> trans,
                                  const Eigen::Ref<const Eigen::Vector4d> quat,
                                  const long& index) const
{
  BoxAbovePlan::diffTrans(res, trans, quat, d_, normal_, index);
}

void BoxAboveFixedPlan::diffQuat(Eigen::Ref<Eigen::Matrix<double, 1, 4>> res,
                                 const Eigen::Ref<const Eigen::Vector3d> trans,
                                 const Eigen::Ref<const Eigen::Vector4d> quat,
                                 const long& index) const
{
  BoxAbovePlan::diffQuat(res, trans, quat, d_, normal_, index);
}

void BoxAboveFixedPlan::fillLinCstr(
    Eigen::Ref<Eigen::Matrix<double, 8, 1>> lb,
    Eigen::Ref<Eigen::Matrix<double, 8, 3>> C) const
{
  lb << -box().vertexMat().transpose() * normal_;
  for (long i = 0; i < 8; ++i) C.row(i) << normal_.transpose();
}

} /* feettrajectory */
