#include <iostream>
#include <limits>

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

void BoxAboveFixedPlan::fillLinCstr(double& lb,
                                    Eigen::Ref<Eigen::MatrixXd> C) const
{
  // Constraint is:
  // for all i in [1:8]
  // -vi.n + d <= n.b
  // Simplified into
  // max_i -vi.n + d <= n.b
  double maxLB = -std::numeric_limits<double>::infinity();
  double val;
  for (size_t i = 0; i < box().vertex().size(); i++)
  {
    val = -box().vertex(i).dot(normal_) + d_;
    if (val > maxLB) maxLB = val;
  }
  lb = maxLB;
  C << normal_.transpose();
}

} /* feettrajectory */
