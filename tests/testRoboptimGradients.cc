#include <iostream>
#include <cmath>

#include <boost/shared_ptr.hpp>

#include <roboptim/core/twice-differentiable-function.hh>
#include <roboptim/core/io.hh>
#include <roboptim/core/numeric-linear-function.hh>
#include <roboptim/core/linear-function.hh>
#include <roboptim/core/optimization-logger.hh>
#include <roboptim/core/solver.hh>
#include <roboptim/core/solver-factory.hh>
#include <roboptim/core/finite-difference-gradient.hh>

#include <manifolds/CartesianProduct.h>
#include <manifolds/Point.h>

#include <cube-stacks/RoboptimCubeStackProblem.hh>
#include <cube-stacks/CubeStackProblemOnR.hh>
#include <cube-stacks/utils/ProblemConfig.hh>
#include <cube-stacks/functions/TotalCost.hh>
#include <cube-stacks/functions/CubeAboveFixedPlanRoboptim.hh>
#include <cube-stacks/functions/CubeAbovePlanRoboptim.hh>
#include <cube-stacks/functions/Norm1Quaternion.hh>
#include <cube-stacks/functions/Norm1Vector.hh>

using namespace cubestacks;
using namespace roboptim;

typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericFunction<T> Function_t;
typedef roboptim::GenericLinearFunction<T> LinearFunction_t;
typedef roboptim::GenericDifferentiableFunction<T> DiffFunction_t;
typedef roboptim::Solver<T> solver_t;
typedef roboptim::OptimizationLogger<solver_t> logger_t;

void displayGradient(
    const DiffFunction_t& function,
    const typename DiffFunction_t::vector_t& x,
    typename DiffFunction_t::size_type i)
{
  GenericFiniteDifferenceGradient<T> fdfunction(function);
  typename GenericFiniteDifferenceGradient<T>::gradient_t grad =
      function.gradient(x, i);
  typename GenericFiniteDifferenceGradient<T>::gradient_t fdgrad =
      fdfunction.gradient(x, i);

  std::cout << "#" << grad << std::endl
            << "#" << fdgrad << std::endl;
}

int main(void)
{
  std::srand(static_cast<uint>(std::time(0)));
  std::string testDir = TESTS_DATA_DIR;
  std::string ymlPath = testDir + "/stackCubesRoboptim.yml";
  ProblemConfig config(ymlPath);


  int nCubes = config["nCubes"];
  IndexManager indexManager(nCubes);

  mnf::CartesianProduct* M = CubeStackProblemOnR::buildManifold(nCubes);

  M->display();

  RoboptimCubeStackProblem myProb(*M, ymlPath);
  Eigen::VectorXd v0 = myProb.pgsProb().findInitPoint();
  
  std::cout << "Norm1Vectors:" << std::endl;
  for (size_t i = 0; i < myProb.Norm1Vectors().size(); i++) 
  {
    displayGradient(*(myProb.Norm1Vectors()[i]), v0, 0);
  }
  std::cout << "Norm1Quaternions:" << std::endl;
  for (size_t i = 0; i < myProb.Norm1Quaternions().size(); i++) 
  {
    displayGradient(*(myProb.Norm1Quaternions()[i]), v0, 0);
  }
  std::cout << "cubeAboveFixedPlanCstrs:" << std::endl;
  for (size_t i = 0; i < myProb.cubeAboveFixedPlanCstrs().size(); i++) 
  {
    displayGradient(*(myProb.cubeAboveFixedPlanCstrs()[i]), v0, 0);
  }
  std::cout << "cubeAbovePlanCstrs:" << std::endl;
  for (size_t i = 0; i < myProb.cubeAbovePlanCstrs().size(); i++) 
  {
    displayGradient(*(myProb.cubeAbovePlanCstrs()[i]), v0, 0);
  }

  return 0;
}
