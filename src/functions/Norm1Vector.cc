#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cube-stacks/functions/Norm1Vector.hh>

namespace cubestacks
{
using namespace Eigen;

Norm1Vector::Norm1Vector(const IndexManager& indexMngr, const int& planIndex)
    : DiffFunction_t(indexMngr.totalDim(), 1, "Norm-1 Vector"),
      planIndex_(planIndex),
      indexManager_(indexMngr)
{
}

Norm1Vector::~Norm1Vector() {}

void Norm1Vector::impl_compute(result_ref res, const_argument_ref x) const
{
  Eigen::Vector3d n = indexManager_.getPlanNormal(planIndex_, x);
  res[0] = n.squaredNorm() - 1;
}

void Norm1Vector::impl_gradient(gradient_ref grad, const_argument_ref x,
                                size_type) const
{
  Eigen::Vector3d n = indexManager_.getPlanNormal(planIndex_, x);
  int nBegin = indexManager_.getPlaneNormalBegin(planIndex_);
  grad[nBegin + 0] = 2 * n[0];
  grad[nBegin + 1] = 2 * n[1];
  grad[nBegin + 2] = 2 * n[2];
}
}  // end of namespace inertialidentification
