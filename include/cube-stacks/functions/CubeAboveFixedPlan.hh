#pragma once

#include <iostream>
#include <Eigen/Core>

#include <cube-stacks/utils/Cube.hh>
#include <cube-stacks/functions/CubeAbovePlan.hh>

namespace cubestacks
{
class CubeAboveFixedPlan : public CubeAbovePlan
{
 public:
  CubeAboveFixedPlan(const Cube& c, const Eigen::Vector3d& normal,
                     const double& d);
  virtual ~CubeAboveFixedPlan();

  void compute(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res,
               const Eigen::Ref<const Eigen::Vector3d> trans,
               const Eigen::Ref<const Eigen::Vector4d> quat) const;
  void diffTrans(Eigen::Ref<Eigen::Matrix<double, 8, 3>> res,
                 const Eigen::Ref<const Eigen::Vector3d> trans,
                 const Eigen::Ref<const Eigen::Vector4d> quat) const;
  void diffQuat(Eigen::Ref<Eigen::Matrix<double, 8, 4>> res,
                const Eigen::Ref<const Eigen::Vector3d> trans,
                const Eigen::Ref<const Eigen::Vector4d> quat) const;

  void diffTrans(Eigen::Ref<Eigen::Matrix<double, 1, 3>> res,
                 const Eigen::Ref<const Eigen::Vector3d> trans,
                 const Eigen::Ref<const Eigen::Vector4d> quat,
                 const long& index) const;
  void diffQuat(Eigen::Ref<Eigen::Matrix<double, 1, 4>> res,
                const Eigen::Ref<const Eigen::Vector3d> trans,
                const Eigen::Ref<const Eigen::Vector4d> quat,
                const long& index) const;

 private:
  Eigen::Vector3d normal_;
  double d_;
};

}  // end of namespace cubestacks

