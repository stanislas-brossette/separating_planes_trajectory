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

#include <cube-stacks/functions/Cost.hh>
#include <cube-stacks/functions/CubeAboveFixedPlan.hh>
#include <cube-stacks/utils/Cube.hh>
#include <cube-stacks/utils/Plan.hh>

#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#define DESC_CUBESTACK_MANIFOLD(z, n, data)\
  BOOST_PP_COMMA_IF(n)ROBOPTIM_REAL_SPACE(3), roboptim::SO3Quat

#define DESC_CUBESTACK_PLAN_MANIFOLD(z, n, data)\
  , roboptim::S2, ROBOPTIM_REAL_SPACE(1)

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

//class CubeStackProblem3Cube
//{
  //ROBOPTIM_DESC_MANIFOLD(R3xSO3, ROBOPTIM_REAL_SPACE(3), roboptim::SO3Quat);
  //ROBOPTIM_DESC_MANIFOLD(PROBLEM_MANIFOLD, BOOST_PP_REPEAT(3, DESC_CUBESTACK_MANIFOLD, ""));
//};

//#pragma comment "LOL"
//#define CUBESTACKPROBLEM(N)\
//struct CubeStackProblemCube##N\
//{\
//ROBOPTIM_DESC_MANIFOLD(R3xSO3, ROBOPTIM_REAL_SPACE(3), roboptim::SO3Quat);\
//ROBOPTIM_DESC_MANIFOLD(PROBLEM_MANIFOLD, BOOST_PP_REPEAT(N, DESC_CUBESTACK_MANIFOLD, "") BOOST_PP_REPEAT(N*(N-1)/2, DESC_CUBESTACK_PLAN_MANIFOLD, ""));\
//};
//
//template<typename T>
//class CubeStackProblemNCube
//{
//  ROBOPTIM_NAMED_FUNCTION_BINDING(Cost_On_M, Cost, T);
//};
//
//CUBESTACKPROBLEM(5)
//auto pb = CubeStackProblemNCube<CubeStackProblemCube5::PROBLEM_MANIFOLD>();

class CubeStackProblem1Cube
{
  ROBOPTIM_DESC_MANIFOLD(R3, ROBOPTIM_REAL_SPACE(3));
  ROBOPTIM_DESC_MANIFOLD(SO3, roboptim::SO3Quat);
  ROBOPTIM_DESC_MANIFOLD(R3xSO3, ROBOPTIM_REAL_SPACE(3), roboptim::SO3Quat);
  //ROBOPTIM_DESC_MANIFOLD(R3xSO3, R3, SO3);

  ROBOPTIM_NAMED_FUNCTION_BINDING(Cost_On_M, Cost, R3xSO3);
  ROBOPTIM_NAMED_FUNCTION_BINDING(CubeAboveF_On_M, CubeAboveFixedPlan, R3xSO3);
public:
  CubeStackProblem1Cube ();
  virtual ~CubeStackProblem1Cube ();
  void init(const Eigen::VectorXd& v, const double& ramdomCoeff = 0);
  //void initTangent(const Eigen::VectorXd& v);
  void solve(const std::string& solvername);
  const mnf::CartesianProduct& M();

  void logdir(const std::string& logdir);

  size_t nCubes_;
  size_t nPlanes_;
  double wallZPlus_;
  double wallXPlus_;
  double wallXMinus_;
  double wallYPlus_;
  double wallYMinus_;
  Eigen::Vector3d normalZPlus_;
  Eigen::Vector3d normalXPlus_;
  Eigen::Vector3d normalXMinus_;
  Eigen::Vector3d normalYPlus_;
  Eigen::Vector3d normalYMinus_;

private:
  std::unique_ptr<roboptim::ProblemOnManifold<T>> manifoldProblem_;
  std::vector<Cube> cubes_;
  std::vector<Plan> plans_;

  std::shared_ptr<Cost_On_M> cost_;
  std::vector<std::shared_ptr<CubeAboveF_On_M>> cstrFixPlan_;

  mnf::RealSpace R3_;
  mnf::SO3<mnf::ExpMapQuaternion> SO3_;
  mnf::CartesianProduct M_;

  roboptim::ManifoldProblemFactory<T> probFactory_;
  //ProblemConfig config_;
  std::string logdir_;
};
} /* cubestacks */
