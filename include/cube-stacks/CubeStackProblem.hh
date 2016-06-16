#pragma once

#include <iostream>

#include <Eigen/Core>

#include <roboptim/core/linear-function.hh>
#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/optimization-logger.hh>
#include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>
#include <roboptim/core/manifold-map/decorator/manifold-problem-factory.hh>

#include <manifolds/RealSpace.h>
#include <manifolds/SO3.h>
#include <manifolds/S2.h>
#include <manifolds/CartesianProduct.h>

namespace cubestacks
{
using namespace Eigen;

typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericFunction<T> Function_t;
typedef roboptim::GenericLinearFunction<T> LinearFunction_t;
typedef roboptim::Solver<T> solver_t;
typedef roboptim::OptimizationLogger<solver_t> logger_t;

typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

class CubeStackProblem
{
  ROBOPTIM_DESC_MANIFOLD(R1, ROBOPTIM_REAL_SPACE(1));
  ROBOPTIM_DESC_MANIFOLD(R3, ROBOPTIM_REAL_SPACE(3));
  ROBOPTIM_DESC_MANIFOLD(SO3, roboptim::SO3);
  ROBOPTIM_DESC_MANIFOLD(S2, roboptim::S2);
  ROBOPTIM_DESC_MANIFOLD(R3xSO3, ROBOPTIM_REAL_SPACE(3), roboptim::SO3);
public:
  CubeStackProblem ();
  virtual ~CubeStackProblem ();
  std::unique_ptr<roboptim::ProblemOnManifold<T>> manifoldProblem_;

private:
  //std::vector<Cube>
};
} /* cubestacks */ 
