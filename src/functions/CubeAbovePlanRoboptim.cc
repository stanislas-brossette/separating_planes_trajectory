#include <iostream>

#include <cube-stacks/functions/CubeAbovePlanRoboptim.hh>
#include <cube-stacks/utils/IndexManager.hh>

namespace cubestacks
{
CubeAbovePlanRoboptim::CubeAbovePlanRoboptim(const IndexManager& indexManager,
                                             CubeAbovePlan& c,
                                             const int& planIndex,
                                             const bool& above)
    : DiffFunction_t(indexManager.totalDim(), 8, "CubeAbovePlanRoboptim"),
      c_(c),
      cubeIndex_(c_.cube().index()),
      planIndex_(planIndex),
      cubeAbove_(above),
      coeff_(cubeAbove_ ? 1 : -1),
      indexManager_(indexManager)
{
}

CubeAbovePlanRoboptim::~CubeAbovePlanRoboptim() {}

void CubeAbovePlanRoboptim::impl_compute(result_ref res,
                                         const_argument_ref x) const
{
  // std::cout << "CubeAbovePlanRoboptim::impl_compute" << std::endl;
  // std::cout << "x: " << x.transpose() << std::endl;
  Eigen::Vector3d trans = indexManager_.getCubeTrans(cubeIndex_, x);
  Eigen::Vector4d quat = indexManager_.getCubeQuat(cubeIndex_, x);
  double d = coeff_ * indexManager_.getPlanDist(planIndex_, x);
  Eigen::Vector3d n = coeff_ * indexManager_.getPlanNormal(planIndex_, x);
  c_.compute(res, trans, quat, d, n);
  // std::cout << "resCubes: " << res.transpose() << std::endl;
}

void CubeAbovePlanRoboptim::impl_gradient(gradient_ref grad,
                                          const_argument_ref x,
                                          size_type index) const
{
  Eigen::Vector3d trans = indexManager_.getCubeTrans(cubeIndex_, x);
  Eigen::Vector4d quat = indexManager_.getCubeQuat(cubeIndex_, x);
  double d = coeff_ * indexManager_.getPlanDist(planIndex_, x);
  Eigen::Vector3d n = coeff_ * indexManager_.getPlanNormal(planIndex_, x);

  Eigen::Matrix<double, 1, 3> diffT;
  Eigen::Matrix<double, 1, 4> diffQ;
  Eigen::Matrix<double, 1, 1> diffD;
  Eigen::Matrix<double, 1, 3> diffN;

  c_.diffTrans(diffT, trans, quat, d, n, index);
  c_.diffQuat(diffQ, trans, quat, d, n, index);
  c_.diffD(diffD, trans, quat, d, n, index);
  c_.diffNormal(diffN, trans, quat, d, n, index);

  int cubeQuatBegin = indexManager_.getCubeQuatBegin(cubeIndex_);
  int cubeTransBegin = indexManager_.getCubeTransBegin(cubeIndex_);
  int planeDistBegin = indexManager_.getPlaneDistBegin(planIndex_);
  int planeNormBegin = indexManager_.getPlaneNormalBegin(planIndex_);

  grad[cubeTransBegin + 0] = diffT(0, 0);
  grad[cubeTransBegin + 1] = diffT(0, 1);
  grad[cubeTransBegin + 2] = diffT(0, 2);

  grad[cubeQuatBegin + 0] = diffQ(0, 0);
  grad[cubeQuatBegin + 1] = diffQ(0, 1);
  grad[cubeQuatBegin + 2] = diffQ(0, 2);
  grad[cubeQuatBegin + 3] = diffQ(0, 3);

  grad[planeDistBegin + 0] = coeff_ * diffD(0, 0);

  grad[planeNormBegin + 0] = coeff_ * diffN(0, 0);
  grad[planeNormBegin + 1] = coeff_ * diffN(0, 1);
  grad[planeNormBegin + 2] = coeff_ * diffN(0, 2);
  // std::cout << "grad: " << grad << std::endl;
}
} /* cubestacks */
