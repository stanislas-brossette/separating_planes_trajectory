#pragma once

#include <Eigen/Core>

namespace cubestacks
{
Eigen::Matrix3d quat2mat(const double& qw, const double& qx, const double& qy,
                         const double& qz);
Eigen::Matrix3d quat2mat(const Eigen::Ref<const Eigen::Vector4d> q);
} /* cubestacks */
