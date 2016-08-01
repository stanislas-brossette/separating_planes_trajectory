#include <iostream>
#include <fstream>
#include <limits>

#include <Eigen/Geometry>

#include <roboptim/core.hh>

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

namespace cubestacks
{

typedef roboptim::EigenMatrixDense T;
typedef roboptim::GenericFunction<T> Function_t;
typedef roboptim::GenericLinearFunction<T> LinearFunction_t;
typedef roboptim::Solver<T> solver_t;
typedef roboptim::OptimizationLogger<solver_t> logger_t;

class solutionVisitor : public boost::static_visitor<Eigen::VectorXd>
{
 public:
  Eigen::VectorXd operator()(roboptim::Result res) const
  {
    std::cout << "Success" << std::endl;
    return res.x;
  }
  //Eigen::VectorXd operator()(roboptim::ResultWithWarnings res) const
  //{
    //std::cout << "Warning" << std::endl;
    //return res.x;
  //}
  Eigen::VectorXd operator()(roboptim::NoSolution) const
  {
    std::cout << "NoSolution" << std::endl;
    Eigen::VectorXd res;
    return res;
  }
  Eigen::VectorXd operator()(roboptim::SolverError) const
  {
    std::cout << "SolverError" << std::endl;
    Eigen::VectorXd res;
    return res;
  }
};

RoboptimCubeStackProblem::RoboptimCubeStackProblem(
    const mnf::Manifold& M, const std::string& configPath)
    : config_(configPath), myProb_(M, configPath)
{
  Eigen::VectorXd v0(M.representationDim());
  v0 = myProb_.findInitPoint();
  M.createRandomPoint(v0);

  int nCubes = config_["nCubes"];
  IndexManager indexManager(nCubes);
  f_ = boost::make_shared<TotalCost>(indexManager);
  solver_t::problem_t problem(f_);

  for (size_t i = 0; i < myProb_.cubeAboveFixedPlanCstrs().size(); i++)
  {
    boost::shared_ptr<CubeAboveFixedPlanRoboptim> g =
        boost::make_shared<CubeAboveFixedPlanRoboptim>(
            indexManager, myProb_.cubeAboveFixedPlanCstrs()[i]);
    cubeAboveFixedPlanCstrs_.push_back(g);
    std::vector<double> scalings;
    typename CubeAboveFixedPlanRoboptim::intervals_t bounds;
    for (int i = 0; i < 8; ++i)
    {
      bounds.push_back(roboptim::Function::makeLowerInterval(0));
      scalings.push_back(1.0);
    }
    problem.addConstraint(g, bounds, scalings);
  }

  for (int iPlan = 0; iPlan < static_cast<int>(myProb_.plans().size()); iPlan++)
  {
    auto iCubeAbove = myProb_.plans()[static_cast<size_t>(iPlan)].cubeAbove();
    auto iCubeBelow = myProb_.plans()[static_cast<size_t>(iPlan)].cubeBelow();

    boost::shared_ptr<CubeAbovePlanRoboptim> gAbove =
        boost::make_shared<CubeAbovePlanRoboptim>(
            indexManager, myProb_.cubeAbovePlanCstrs()[iCubeAbove], iPlan, true);
    cubeAbovePlanCstrs_.push_back(gAbove);
    std::vector<double> scalings;
    typename CubeAboveFixedPlanRoboptim::intervals_t bounds;
    for (int i = 0; i < 8; ++i)
    {
      bounds.push_back(roboptim::Function::makeLowerInterval(0));
      scalings.push_back(1.0);
    }
    problem.addConstraint(gAbove, bounds, scalings);

    boost::shared_ptr<CubeAbovePlanRoboptim> gBelow =
        boost::make_shared<CubeAbovePlanRoboptim>(
            indexManager, myProb_.cubeAbovePlanCstrs()[iCubeBelow], iPlan,
            false);
    cubeAbovePlanCstrs_.push_back(gBelow);
    scalings.clear();
    bounds.clear();
    for (int i = 0; i < 8; ++i)
    {
      bounds.push_back(roboptim::Function::makeLowerInterval(0));
      scalings.push_back(1.0);
    }
    problem.addConstraint(gBelow, bounds, scalings);
  }

  for (int i = 0; i < nCubes; i++)
  {
    boost::shared_ptr<Norm1Quaternion> g =
        boost::make_shared<Norm1Quaternion>(indexManager, i);
    Norm1Quaternions_.push_back(g);
    problem.addConstraint(g, roboptim::Function::makeInterval(0, 0), 1.0);
  }

  for (int i = 0; i < indexManager.nPlanes(); i++)
  {
    boost::shared_ptr<Norm1Vector> g =
        boost::make_shared<Norm1Vector>(indexManager, i);
    Norm1Vectors_.push_back(g);
    problem.addConstraint(g, roboptim::Function::makeInterval(0, 0), 1.0);
  }

  std::cout << problem << std::endl;

  /**********
   *  Init  *
   **********/
  TotalCost::argument_t x (indexManager.totalDim());
  x << v0;
  std::cout << "x: " << x << std::endl;
  problem.startingPoint () = x;
  std::cout << "problem.startingPoint(): " << problem.startingPoint() << std::endl;
  std::string solverName(config_["solverName"]);
  std::cout << "solverName: " << solverName << std::endl;
  roboptim::SolverFactory<solver_t> solverFactory(solverName,
                                                  problem);
  solver_t& solver = solverFactory();

  /************
   *  Result  *
   ************/
  solver_t::result_t res(solver.minimum());
  std::cout << "res: " << res << std::endl;
  Eigen::VectorXd xSol = boost::apply_visitor(solutionVisitor(), res);

  std::cout << "xSol = \n" << xSol << std::endl;
  mnf::Point pSol = M.createPoint(xSol);
  std::cout << "pSol: " << pSol << std::endl;
  myProb_.fileForMatlab(
      "/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      "stackCubesCFSQP.m",
      pSol);
}

RoboptimCubeStackProblem::~RoboptimCubeStackProblem() {}



}
