#include <cube-stacks/CubeStackProblem1Cube.hh>

#include <roboptim/core/solver-factory.hh>
#include <roboptim/core/solver.hh>

namespace cubestacks
{
CubeStackProblem1Cube::CubeStackProblem1Cube()
  : R3_(3),
  M_({&R3_, &SO3_})
{
  cubes_.push_back(Cube(1.1));

  std::cout << "M_: " << M_ << std::endl;
  /*******************
   *  Cost Function  *
   *******************/
  cost_ = std::make_shared<Cost_On_M>();
  probFactory_.addObjective(cost_, M_);

  /*****************
   *  Constraints  *
   *****************/
  cstrFixPlan_.push_back(std::make_shared<CubeAboveF_On_M>(cubes_[0], normalZPlus_, wallZPlus_));
  cstrFixPlan_.push_back(std::make_shared<CubeAboveF_On_M>(cubes_[0], normalXMinus_, wallXMinus_));
  cstrFixPlan_.push_back(std::make_shared<CubeAboveF_On_M>(cubes_[0], normalXPlus_, wallXPlus_));
  cstrFixPlan_.push_back(std::make_shared<CubeAboveF_On_M>(cubes_[0], normalYMinus_, wallYMinus_));
  cstrFixPlan_.push_back(std::make_shared<CubeAboveF_On_M>(cubes_[0], normalYPlus_, wallYPlus_));

  std::vector<double> scalings;
  typename CubeAboveFixedPlan::intervals_t boundsNorm;
  for (int i = 0; i < 8; ++i)
  {
    boundsNorm.push_back(roboptim::Function::makeLowerInterval(0));
    scalings.push_back(1.0);
  }
  for (size_t i = 0; i < cstrFixPlan_.size(); i++) 
  {
    probFactory_.addConstraint(cstrFixPlan_[i], M_).setBounds(boundsNorm).setScaling(scalings);
  }
  scalings.clear();
  /*************
   *   Bounds  *
   *************/

  /***************************
   *  Problem Instanciation  *
   ***************************/
  manifoldProblem_ = probFactory_.getProblem();
}

CubeStackProblem1Cube::~CubeStackProblem1Cube()
{
}

void CubeStackProblem1Cube::logdir(const std::string& logdir)
{
  logdir_ = logdir;
}

void CubeStackProblem1Cube::init(const Eigen::VectorXd& v, const double& randomCoeff)
{
  /**********
   *  Init  *
   **********/
  Eigen::VectorXd w = v;
  manifoldProblem_->getManifold().forceOnM(w, v);
  mnf::Point x0 = manifoldProblem_->getManifold().createPoint(w);
  if (randomCoeff != 0)
    x0.increment(randomCoeff * Eigen::VectorXd::Random(manifoldProblem_->getManifold().tangentDim()));
  manifoldProblem_->startingPoint() = x0.value();
}

void CubeStackProblem1Cube::solve(const std::string& solverName)
{
  /************
   *  Solver  *
   ************/
  roboptim::SolverFactory<solver_t> solverFactory(solverName,
                                                  *manifoldProblem_);
  solver_t& solver = solverFactory();

  //std::shared_ptr<logger_t> logger;
  //if(config_["log_roboptim"])
    //logger = std::make_shared<logger_t>(solver, logdir_);

  //setAllPGSolverOptions(solver, config_);
  std::cout << *manifoldProblem_ << std::endl;
  //std::cout << solver << std::endl;

  /************
   *  Result  *
   ************/
  solver_t::result_t res(solver.minimum());
  std::cout << res << std::endl;
}

const mnf::CartesianProduct& CubeStackProblem1Cube::M()
{
  return M_;
}

} /* cubestack */ 
