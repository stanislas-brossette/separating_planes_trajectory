#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cube-stacks/functions/Norm1Quaternion.hh>

namespace cubestacks
{
using namespace Eigen;

Norm1Quaternion::Norm1Quaternion(const IndexManager& indexMngr,
                                 const int& cubeIndex)
    : DiffFunction_t(indexMngr.totalDim(), 1, "Norm-1 quaternion"),
      cubeIndex_(cubeIndex),
      indexManager_(indexMngr)
{
}

Norm1Quaternion::~Norm1Quaternion() {}

void Norm1Quaternion::impl_compute(result_ref res, const_argument_ref x) const
{
  //std::cout << "\nNorm1Quaternion::impl_compute" << std::endl;
  //std::cout << "x: " << x.transpose() << std::endl;
  Eigen::Vector4d quat = indexManager_.getCubeQuat(cubeIndex_, x);
  res[0] = quat.squaredNorm() - 1;
  //std::cout << "resQuat: " << res[0] << std::endl;
}

void Norm1Quaternion::impl_gradient(gradient_ref grad, const_argument_ref x,
                                    size_type) const
{
  //std::cout << "Norm1Quaternion::impl_gradient" << std::endl;
  Eigen::Vector4d q = indexManager_.getCubeQuat(cubeIndex_, x);
  int qBegin = indexManager_.getCubeQuatBegin(cubeIndex_);
  grad[qBegin + 0] = 2 * q[0];
  grad[qBegin + 1] = 2 * q[1];
  grad[qBegin + 2] = 2 * q[2];
  grad[qBegin + 3] = 2 * q[3];
  //std::cout << "grad: " << grad << std::endl;
}
}  // end of namespace inertialidentification
