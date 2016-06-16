#pragma once

#include <Eigen/Core>

#include <roboptim/core/differentiable-function.hh>

namespace cubestacks
{
typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericDifferentiableFunction<T> DiffFunction_t;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

/// @brief Class representing the cost function
/// cost = sum(z_i^2)

class Cost : public DiffFunction_t
{
public:
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_(
      roboptim::GenericDifferentiableFunction<T>);

  Cost();

  // Considers that the input comes in the following order:
  // x = [x0, y0, z0, x1, ..., zn] 
  void impl_compute(result_ref res, const_argument_ref x) const;

  void impl_gradient(gradient_ref grad, const_argument_ref, size_type ) const;

  virtual ~Cost();

private:
};

}  // end of namespace cubestacks

