#include <iostream>

#include <cube-stacks/functions/CubeAboveFixedPlanRoboptim.hh>
#include <cube-stacks/utils/IndexManager.hh>

namespace cubestacks
{
CubeAboveFixedPlanRoboptim::CubeAboveFixedPlanRoboptim(const IndexManager& indexManager, CubeAboveFixedPlan& c)
    : DiffFunction_t(indexManager.totalDim(), 8, "CubeAboveFixedPlanRoboptim"),
      c_(c),
      cubeIndex_(c_.cube().index()),
      indexManager_(indexManager)
{
}

CubeAboveFixedPlanRoboptim::~CubeAboveFixedPlanRoboptim() {}

void CubeAboveFixedPlanRoboptim::impl_compute(result_ref res, const_argument_ref x) const
{
  //std::cout << "CubeAboveFixedPlanRoboptim::impl_compute" << std::endl;
  //std::cout << "x: " << x.transpose() << std::endl;
  Eigen::Vector3d trans = indexManager_.getCubeTrans(cubeIndex_, x);
  Eigen::Vector4d quat = indexManager_.getCubeQuat(cubeIndex_, x);
  c_.compute(res, trans, quat);
  //std::cout << "resCubes: " << res.transpose() << std::endl;
}

void CubeAboveFixedPlanRoboptim::impl_gradient(gradient_ref grad,
                                                       const_argument_ref x,
                                                       size_type index) const
{
  //std::cout << "CubeAboveFixedPlanRoboptim::impl_gradient" << std::endl;
  Eigen::Vector3d trans = indexManager_.getCubeTrans(cubeIndex_, x);
  Eigen::Vector4d quat = indexManager_.getCubeQuat(cubeIndex_, x);
  Eigen::Matrix<double, 8, 3> diffT;
  Eigen::Matrix<double, 8, 4> diffQ;
  c_.diffTrans(diffT, trans, quat);
  c_.diffQuat(diffQ, trans, quat);
  int cubeQuatBegin = indexManager_.getCubeQuatBegin(cubeIndex_);
  int cubeTransBegin = indexManager_.getCubeTransBegin(cubeIndex_);
  grad[cubeTransBegin + 0] = diffT(index,0);
  grad[cubeTransBegin + 1] = diffT(index,1);
  grad[cubeTransBegin + 2] = diffT(index,2);
  grad[cubeQuatBegin + 0] = diffQ(index,0);
  grad[cubeQuatBegin + 1] = diffQ(index,1);
  grad[cubeQuatBegin + 2] = diffQ(index,2);
  grad[cubeQuatBegin + 3] = diffQ(index,3);
  //std::cout << "grad: " << grad << std::endl;
}
} /* cubestacks */
