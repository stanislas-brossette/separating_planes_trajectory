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

typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericFunction<T> Function_t;
typedef roboptim::GenericLinearFunction<T> LinearFunction_t;
typedef roboptim::Solver<T> solver_t;
typedef roboptim::OptimizationLogger<solver_t> logger_t;

//class solutionVisitor : public boost::static_visitor<Eigen::VectorXd>
//{
 //public:
  //Eigen::VectorXd operator()(roboptim::Result res) const
  //{
    //std::cout << "Success" << std::endl;
    //return res.x;
  //}
  ////Eigen::VectorXd operator()(roboptim::ResultWithWarnings res) const
  ////{
    ////std::cout << "Warning" << std::endl;
    ////return res.x;
  ////}
  //Eigen::VectorXd operator()(roboptim::NoSolution) const
  //{
    //std::cout << "NoSolution" << std::endl;
    //Eigen::VectorXd res;
    //return res;
  //}
  //Eigen::VectorXd operator()(roboptim::SolverError) const
  //{
    //std::cout << "SolverError" << std::endl;
    //Eigen::VectorXd res;
    //return res;
  //}
//};

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

  return 0;
}
