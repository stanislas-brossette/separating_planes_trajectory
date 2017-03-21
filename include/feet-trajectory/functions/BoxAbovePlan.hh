#pragma once
#include <iostream>
#include <Eigen/Core>
#include <feet-trajectory/utils/Box.hh>

namespace feettrajectory
{
class BoxAbovePlan
{
 public:
  BoxAbovePlan(const Box& c);
  virtual ~BoxAbovePlan();
  void compute(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res,
               const Eigen::Ref<const Eigen::Vector3d> trans,
               const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
               const Eigen::Ref<const Eigen::Vector3d> normal) const;
  void diffTrans(Eigen::Ref<Eigen::Matrix<double, 8, 3>> res,
                 const Eigen::Ref<const Eigen::Vector3d> trans,
                 const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
                 const Eigen::Ref<const Eigen::Vector3d> normal) const;
  void diffQuat(Eigen::Ref<Eigen::Matrix<double, 8, 4>> res,
                const Eigen::Ref<const Eigen::Vector3d> trans,
                const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
                const Eigen::Ref<const Eigen::Vector3d> normal) const;
  void diffD(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res,
             const Eigen::Ref<const Eigen::Vector3d> trans,
             const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
             const Eigen::Ref<const Eigen::Vector3d> normal) const;
  void diffNormal(Eigen::Ref<Eigen::Matrix<double, 8, 3>> res,
                  const Eigen::Ref<const Eigen::Vector3d> trans,
                  const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
                  const Eigen::Ref<const Eigen::Vector3d> normal) const;

  void diffTrans(Eigen::Ref<Eigen::Matrix<double, 1, 3>> res,
                 const Eigen::Ref<const Eigen::Vector3d> trans,
                 const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
                 const Eigen::Ref<const Eigen::Vector3d> normal,
                 const long& index) const;
  void diffQuat(Eigen::Ref<Eigen::Matrix<double, 1, 4>> res,
                const Eigen::Ref<const Eigen::Vector3d> trans,
                const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
                const Eigen::Ref<const Eigen::Vector3d> normal,
                const long& index) const;
  void diffD(Eigen::Ref<Eigen::Matrix<double, 1, 1>> res,
             const Eigen::Ref<const Eigen::Vector3d> trans,
             const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
             const Eigen::Ref<const Eigen::Vector3d> normal,
             const long& index) const;
  void diffNormal(Eigen::Ref<Eigen::Matrix<double, 1, 3>> res,
                  const Eigen::Ref<const Eigen::Vector3d> trans,
                  const Eigen::Ref<const Eigen::Vector4d> quat, const double& d,
                  const Eigen::Ref<const Eigen::Vector3d> normal,
                  const long& index) const;
  void LB(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res) const;
  void UB(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res) const;
  const Box& box() const { return box_; };
  long dim() const { return 8; }

  void fillLinCstrNFixed(double& lb,
                         Eigen::Ref<Eigen::MatrixXd> CboxPart,
                         double& CdPart,
                         Eigen::Vector3d normal,
                         bool isBoxBelow = false) const;

 private:
  Box box_;
};

//cstr = nx*(tx - vx*(2*qy^2 + 2*qz^2 - 1) - vy*(2*qw*qz - 2*qx*qy) + vz*(2*qw*qy + 2*qx*qz)) - d + ny*(ty - vy*(2*qx^2 + 2*qz^2 - 1) + vx*(2*qw*qz + 2*qx*qy) - vz*(2*qw*qx - 2*qy*qz)) + nz*(tz - vz*(2*qx^2 + 2*qy^2 - 1) - vx*(2*qw*qy - 2*qx*qz) + vy*(2*qw*qx + 2*qy*qz))
//diff(cstr, tx) = nx
//diff(cstr, ty) = ny
//diff(cstr, tz) = nz
//diff(cstr, qw) = nz*(2*qw*vy - 4*qx*vz + 2*qz*vx) - ny*(2*qw*vz + 4*qx*vy - 2*qy*vx) + nx*(2*qy*vy + 2*qz*vz)
//diff(cstr, qx) = nx*(2*qw*vz + 2*qx*vy - 4*qy*vx) - nz*(2*qw*vx + 4*qy*vz - 2*qz*vy) + ny*(2*qx*vx + 2*qz*vz)
//diff(cstr, qy) = ny*(2*qw*vx + 2*qy*vz - 4*qz*vy) - nx*(2*qw*vy - 2*qx*vz + 4*qz*vx) + nz*(2*qx*vx + 2*qy*vy)
//diff(cstr, qz) = nz*(2*qx*vy - 2*qy*vx) - ny*(2*qx*vz - 2*qz*vx) + nx*(2*qy*vz - 2*qz*vy)
//diff(cstr, nx) = tx - vx*(2*qy^2 + 2*qz^2 - 1) - vy*(2*qw*qz - 2*qx*qy) + vz*(2*qw*qy + 2*qx*qz)
//diff(cstr, ny) = ty - vy*(2*qx^2 + 2*qz^2 - 1) + vx*(2*qw*qz + 2*qx*qy) - vz*(2*qw*qx - 2*qy*qz)
//diff(cstr, nz) = tz - vz*(2*qx^2 + 2*qy^2 - 1) - vx*(2*qw*qy - 2*qx*qz) + vy*(2*qw*qx + 2*qy*qz)
//diff(cstr, d)  = -1

} /* feettrajectory */
