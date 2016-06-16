#include <cube-stacks/functions/CubeAbovePlan.hh>
#include <cube-stacks/utils/quat2mat.hh>

namespace cubestacks
{
CubeAbovePlan::CubeAbovePlan(const Cube& c) : cube_(c) {}

CubeAbovePlan::~CubeAbovePlan() {}

void CubeAbovePlan::compute(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res,
                            const Eigen::Ref<const Eigen::Vector3d> t,
                            const Eigen::Ref<const Eigen::Vector4d> q,
                            const double& d,
                            const Eigen::Ref<const Eigen::Vector3d> n) const
{
  Eigen::Matrix3d rot = quat2mat(q);
  Eigen::Matrix<double, 3, 8> vertexMatWorld;
  for (size_t i = 0; i < 8; i++)
  {
    Eigen::Vector3d v;
    v = t + rot * cube_.vertex(i);
    res[static_cast<long>(i)] = v.dot(n) - d;
  }
}

void CubeAbovePlan::diffTrans(Eigen::Ref<Eigen::Matrix<double, 8, 3>> res,
                              const Eigen::Ref<const Eigen::Vector3d> ,
                              const Eigen::Ref<const Eigen::Vector4d> ,
                              const double& ,
                              const Eigen::Ref<const Eigen::Vector3d> n) const
{
  for (long i = 0; i < 8; i++)
  {
    res(i, 0) = n.x();
    res(i, 1) = n.y();
    res(i, 2) = n.z();
  }
}

void CubeAbovePlan::diffQuat(Eigen::Ref<Eigen::Matrix<double, 8, 4>> res,
                             const Eigen::Ref<const Eigen::Vector3d>,
                             const Eigen::Ref<const Eigen::Vector4d> q, const double&,
                             const Eigen::Ref<const Eigen::Vector3d> n) const
{
  for (long i = 0; i < 8; i++)
  {
    Eigen::Vector3d v = cube_.vertex(i);
    res( i, 0 ) = n.z() * (2 * q.x() * v.y() - 2 * q.y() * v.x()) -
                  n.y() * (2 * q.x() * v.z() - 2 * q.z() * v.x()) +
                  n.x() * (2 * q.y() * v.z() - 2 * q.z() * v.y());
    res( i, 1 ) =
        n.z() * (2 * q.w() * v.y() - 4 * q.x() * v.z() + 2 * q.z() * v.x()) -
        n.y() * (2 * q.w() * v.z() + 4 * q.x() * v.y() - 2 * q.y() * v.x()) +
        n.x() * (2 * q.y() * v.y() + 2 * q.z() * v.z());
    res( i, 2 ) =
        n.x() * (2 * q.w() * v.z() + 2 * q.x() * v.y() - 4 * q.y() * v.x()) -
        n.z() * (2 * q.w() * v.x() + 4 * q.y() * v.z() - 2 * q.z() * v.y()) +
        n.y() * (2 * q.x() * v.x() + 2 * q.z() * v.z());
    res( i, 3 ) =
        n.y() * (2 * q.w() * v.x() + 2 * q.y() * v.z() - 4 * q.z() * v.y()) -
        n.x() * (2 * q.w() * v.y() - 2 * q.x() * v.z() + 4 * q.z() * v.x()) +
        n.z() * (2 * q.x() * v.x() + 2 * q.y() * v.y());
  }
}

void CubeAbovePlan::diffD(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res,
                          const Eigen::Ref<const Eigen::Vector3d>,
                          const Eigen::Ref<const Eigen::Vector4d>, const double&,
                          const Eigen::Ref<const Eigen::Vector3d>) const
{
  for (long i = 0; i < 8; i++)
  {
    res( i, 0 ) = -1;
  }
}

void CubeAbovePlan::diffNormal(Eigen::Ref<Eigen::Matrix<double, 8, 3>> res,
                               const Eigen::Ref<const Eigen::Vector3d> t,
                               const Eigen::Ref<const Eigen::Vector4d> q,
                               const double&,
                               const Eigen::Ref<const Eigen::Vector3d>) const
{
  for (long i = 0; i < 8; i++)
  {
    Eigen::Vector3d v = cube_.vertex(i);
    res( i, 0 ) = t.x() - v.x() * (2 * q.y() * q.y() + 2 * q.z() * q.z() - 1) -
                  v.y() * (2 * q.w() * q.z() - 2 * q.x() * q.y()) +
                  v.z() * (2 * q.w() * q.y() + 2 * q.x() * q.z());
    res( i, 1 ) = t.y() - v.y() * (2 * q.x() * q.x() + 2 * q.z() * q.z() - 1) +
                  v.x() * (2 * q.w() * q.z() + 2 * q.x() * q.y()) -
                  v.z() * (2 * q.w() * q.x() - 2 * q.y() * q.z());
    res( i, 2 ) = t.z() - v.z() * (2 * q.x() * q.x() + 2 * q.y() * q.y() - 1) -
                  v.x() * (2 * q.w() * q.y() - 2 * q.x() * q.z()) +
                  v.y() * (2 * q.w() * q.x() + 2 * q.y() * q.z());
  }
}

void CubeAbovePlan::LB(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res) const
{
  res << 0, 0, 0, 0, 0, 0, 0, 0;
}

void CubeAbovePlan::UB(Eigen::Ref<Eigen::Matrix<double, 8, 1>> res) const
{
  double inf= std::numeric_limits<double>::infinity();
  res << inf, inf, inf, inf, inf, inf, inf, inf;
}

} /* cubestacks */
