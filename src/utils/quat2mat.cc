#include <cube-stacks/utils/quat2mat.hh>

namespace cubestacks
{
Eigen::Matrix3d quat2mat(const double& qx, const double& qy, const double& qz,
                         const double& qw)
{
  Eigen::Matrix3d res;
  res << 1 - 2 * qy* qy - 2 * qz* qz, 2 * qx* qy - 2 * qz* qw,
      2 * qx* qz + 2 * qy* qw, 2 * qx* qy + 2 * qz* qw,
      1 - 2 * qx* qx - 2 * qz* qz, 2 * qy* qz - 2 * qx* qw,
      2 * qx* qz - 2 * qy* qw, 2 * qy* qz + 2 * qx* qw,
      1 - 2 * qx* qx - 2 * qy* qy;
  return res;
}
Eigen::Matrix3d quat2mat(const Eigen::Ref<const Eigen::Vector4d> q)
{
  return quat2mat(q[0], q[1], q[2], q[3]);
}
} /* cubestacks */
