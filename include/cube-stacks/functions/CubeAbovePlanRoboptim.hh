#pragma once

#include <Eigen/Core>

#include <roboptim/core.hh>

#include <cube-stacks/functions/CubeAbovePlan.hh>
#include <cube-stacks/utils/IndexManager.hh>

namespace cubestacks
{
typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericDifferentiableFunction<T> DiffFunction_t;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

class CubeAbovePlanRoboptim : public DiffFunction_t
{
public:
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_(
      roboptim::GenericDifferentiableFunction<T>);

  CubeAbovePlanRoboptim(const IndexManager& indexManager, CubeAbovePlan& c,
                        const int& planIndex, const bool& cubeAbove = true);

  void impl_compute(result_ref res, const_argument_ref x) const;

  void impl_gradient(gradient_ref grad, const_argument_ref, size_type ) const;

  virtual ~CubeAbovePlanRoboptim();

private:
  CubeAbovePlan& c_;
  int cubeIndex_;
  int planIndex_;
  bool cubeAbove_;
  int coeff_;
  IndexManager indexManager_;
};

}  // end of namespace cubestacks


