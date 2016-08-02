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
#include <cube-stacks/utils/IndexManager.hh>
#include <cube-stacks/functions/TotalCost.hh>
#include <cube-stacks/functions/CubeAboveFixedPlanRoboptim.hh>
#include <cube-stacks/functions/CubeAbovePlanRoboptim.hh>
#include <cube-stacks/functions/Norm1Quaternion.hh>
#include <cube-stacks/functions/Norm1Vector.hh>

namespace cubestacks
{
class solutionVisitor
    : public boost::static_visitor<resultRoboptimCubeStackProblem>
{
 public:
  resultRoboptimCubeStackProblem operator()(roboptim::Result res) const
  {
    std::cout << "Success" << std::endl;
    resultRoboptimCubeStackProblem rws;
    rws.status = 1;
    rws.xSol = res.x;
    rws.obj_star = res.value[0];
    return rws;
  }
  // Eigen::VectorXd operator()(roboptim::ResultWithWarnings res) const
  //{
  // std::cout << "Warning" << std::endl;
  // return res.x;
  //}
  resultRoboptimCubeStackProblem operator()(roboptim::NoSolution) const
  {
    std::cout << "NoSolution" << std::endl;
    resultRoboptimCubeStackProblem rws;
    rws.status = 0;
    return rws;
  }
  resultRoboptimCubeStackProblem operator()(roboptim::SolverError) const
  {
    std::cout << "SolverError" << std::endl;
    resultRoboptimCubeStackProblem rws;
    rws.status = 0;
    return rws;
  }
};

RoboptimCubeStackProblem::RoboptimCubeStackProblem(
    const mnf::Manifold& M, const std::string& configPath)
    : nCubes_(static_cast<size_t>(M(0).dim()) / 7),
      nPlans_(nCubes_ * (nCubes_ - 1) / 2),
      config_(configPath),
      pgsProb_(M, configPath),
      indexManager_(static_cast<int>(nCubes_))
{
  Eigen::VectorXd v0(M.representationDim());

  f_ = boost::make_shared<TotalCost>(indexManager_);
  problem_ = std::make_shared<solver_t::problem_t>(f_);

  for (size_t i = 0; i < pgsProb_.cubeAboveFixedPlanCstrs().size(); i++)
  {
    boost::shared_ptr<CubeAboveFixedPlanRoboptim> g =
        boost::make_shared<CubeAboveFixedPlanRoboptim>(
            indexManager_, pgsProb_.cubeAboveFixedPlanCstrs()[i]);
    cubeAboveFixedPlanCstrs_.push_back(g);
    std::vector<double> scalings;
    typename CubeAboveFixedPlanRoboptim::intervals_t bounds;
    for (int i = 0; i < 8; ++i)
    {
      bounds.push_back(roboptim::Function::makeLowerInterval(0));
      scalings.push_back(1.0);
    }
    problem_->addConstraint(g, bounds, scalings);
  }

  for (int iPlan = 0; iPlan < static_cast<int>(pgsProb_.plans().size());
       iPlan++)
  {
    auto iCubeAbove = pgsProb_.plans()[static_cast<size_t>(iPlan)].cubeAbove();
    auto iCubeBelow = pgsProb_.plans()[static_cast<size_t>(iPlan)].cubeBelow();

    boost::shared_ptr<CubeAbovePlanRoboptim> gAbove =
        boost::make_shared<CubeAbovePlanRoboptim>(
            indexManager_, pgsProb_.cubeAbovePlanCstrs()[iCubeAbove], iPlan,
            true);
    cubeAbovePlanCstrs_.push_back(gAbove);
    std::vector<double> scalings;
    typename CubeAboveFixedPlanRoboptim::intervals_t bounds;
    for (int i = 0; i < 8; ++i)
    {
      bounds.push_back(roboptim::Function::makeLowerInterval(0));
      scalings.push_back(1.0);
    }
    problem_->addConstraint(gAbove, bounds, scalings);

    boost::shared_ptr<CubeAbovePlanRoboptim> gBelow =
        boost::make_shared<CubeAbovePlanRoboptim>(
            indexManager_, pgsProb_.cubeAbovePlanCstrs()[iCubeBelow], iPlan,
            false);
    cubeAbovePlanCstrs_.push_back(gBelow);
    scalings.clear();
    bounds.clear();
    for (int i = 0; i < 8; ++i)
    {
      bounds.push_back(roboptim::Function::makeLowerInterval(0));
      scalings.push_back(1.0);
    }
    problem_->addConstraint(gBelow, bounds, scalings);
  }

  for (int i = 0; i < static_cast<int>(nCubes_); i++)
  {
    boost::shared_ptr<Norm1Quaternion> g =
        boost::make_shared<Norm1Quaternion>(indexManager_, i);
    Norm1Quaternions_.push_back(g);
    problem_->addConstraint(g, roboptim::Function::makeInterval(0, 0), 1.0);
  }

  for (int i = 0; i < indexManager_.nPlanes(); i++)
  {
    boost::shared_ptr<Norm1Vector> g =
        boost::make_shared<Norm1Vector>(indexManager_, i);
    Norm1Vectors_.push_back(g);
    problem_->addConstraint(g, roboptim::Function::makeInterval(0, 0), 1.0);
  }
}

void RoboptimCubeStackProblem::init(const Eigen::VectorXd& x0)
{
  TotalCost::argument_t x(indexManager_.totalDim());
  x << x0;
  problem_->startingPoint() = x;
}

void RoboptimCubeStackProblem::setCFSQPparameterInt(solver_t& solver,
                                                    const std::string& s)
{
  if (config_.has(s))
  {
    roboptim::Parameter param;
    param.value = config_[s].asInt();
    solver.parameters()[s] = param;
  }
}

void RoboptimCubeStackProblem::setCFSQPparameterDouble(solver_t& solver,
                                                       const std::string& s)
{
  if (config_.has(s))
  {
    roboptim::Parameter param;
    param.value = config_[s].asDouble();
    solver.parameters()[s] = param;
  }
}

void RoboptimCubeStackProblem::setCFSQPparameters(solver_t& solver)
{
  setCFSQPparameterInt(solver, "max-iterations");
  setCFSQPparameterInt(solver, "cfsqp.mode");
  setCFSQPparameterInt(solver, "cfsqp.iprint");
  setCFSQPparameterDouble(solver, "cfsqp.bigbnd");
  setCFSQPparameterDouble(solver, "cfsqp.eps");
  setCFSQPparameterDouble(solver, "cfsqp.epseqn");
  setCFSQPparameterDouble(solver, "cfsqp.udelta");
  setCFSQPparameterDouble(solver, "cfsqp.objeps");
  setCFSQPparameterDouble(solver, "cfsqp.objrep");
  setCFSQPparameterDouble(solver, "cfsqp.gLgeps");
  setCFSQPparameterInt(solver, "cfsqp.nstop");
}

resultRoboptimCubeStackProblem RoboptimCubeStackProblem::solve(
    const std::string& solverName)
{
  roboptim::SolverFactory<solver_t> solverFactory(solverName, *problem_);
  solver_t& solver = solverFactory();

  nIterations_ = 0;
  solver.setIterationCallback(
      [this](const solver_t::problem_t&, solver_t::solverState_t&)
      {
        nIterations_++;
      });

  if (solverName.compare("cfsqp") == 0)
  {
    setCFSQPparameters(solver);
  }

  //std::cout << "solver: " << solver << std::endl;
  /************
   *  Result  *
   ************/
  solver_t::result_t res(solver.minimum());
  auto result = boost::apply_visitor(solutionVisitor(), res);
  if (result.status)
  {
    std::cout << "============= " << solverName << " Solution after " << nIterations_
              << " iterations ================" << std::endl;
    std::cerr << "============= " << solverName << " Solution after " << nIterations_
              << " iterations ================" << std::endl;
  }
  else
  {
    std::cout << "============= " << solverName << "  Failed  after " << nIterations_
              << " iterations ================" << std::endl;
    std::cerr << "============= " << solverName << "  Failed  after " << nIterations_
              << " iterations ================" << std::endl;
  }

  return result;

  // std::cout << "xSol = \n" << rws.xSol << std::endl;
  // mnf::Point pSol = pgsProb_.M().createPoint(rws.xSol);
  // std::cout << "pSol: " << pSol << std::endl;
  // pgsProb_.fileForMatlab(
  //"/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
  //"stackCubesCFSQP.m",
  // pSol);
}

RoboptimCubeStackProblem::~RoboptimCubeStackProblem() {}

void RoboptimCubeStackProblem::name(const std::string& n) { name_ = n; }
}
