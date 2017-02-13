#include <cstdlib>
#include <ctime>
#include <iostream>

#include <manifolds/CartesianProduct.h>

#include <pgsolver/utils/finiteDiff.h>
#include <pgsolver/utils/loadFromYaml.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverTrustRegionFilter.h>

#include <feet-trajectory/BoxesHullTrajProblem.hh>
#include <feet-trajectory/utils/Printer.hh>

using namespace feettrajectory;

int main(int argc, char *argv[])
{
  std::cout << "Generate Trajectory: Loads a config and generates a trajectory "
               "between initPos and finalPos that avoids the obstacles and "
               "some fixed planes" << std::endl;

  std::srand(static_cast<uint>(std::time(0)));
  std::string ymlPath = std::string(CONFIGS_DATA_DIR) + "/" + argv[1];
  ProblemConfig config(ymlPath);

  int nBoxes = config["nBoxes"];
  int nObstacles = config["obstacles"].asVecBox().size();
  int nFixedPlanes = config["nFixedPlanes"];

  mnf::CartesianProduct* M =
      BoxesHullTrajProblem::buildManifold(nBoxes, nObstacles);
  M->display();

  BoxesHullTrajProblem myProb(*M, ymlPath);

  Eigen::VectorXd v0(M->representationDim());
  v0 = myProb.findInitPoint();
  myProb.M().forceOnM(v0, v0);
  Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[",
                           "]");

  std::cout << "v0: " << v0.format(HeavyFmt) << std::endl;

  pgs::utils::finiteDiffCheck(myProb);

  myProb.M().forceOnM(v0, v0);
  mnf::Point x0 = myProb.M().createPoint(v0);

  // Solver mySolver;
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  pgs::Results res = mySolver.solve();

  std::cout << "xSol = \n" << res.x_star << std::endl;

  print(config["logName"], myProb, res.x_star);

  return 0;
}
