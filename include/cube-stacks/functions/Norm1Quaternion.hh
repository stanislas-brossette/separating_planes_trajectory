#pragma once

#include <Eigen/Core>

#include <roboptim/core.hh>

#include <cube-stacks/utils/IndexManager.hh>

namespace cubestacks
{
typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericDifferentiableFunction<T> DiffFunction_t;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

class Norm1Quaternion : public DiffFunction_t
{
public:
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_(
      roboptim::GenericDifferentiableFunction<T>);

  Norm1Quaternion(const IndexManager& indexMngr, const int& cubeIndex);

  void impl_compute(result_ref res, const_argument_ref x) const;

  void impl_gradient(gradient_ref grad, const_argument_ref, size_type ) const;

  virtual ~Norm1Quaternion();

private:
  int cubeIndex_;
  IndexManager indexManager_;
};

}  // end of namespace cubestacks



