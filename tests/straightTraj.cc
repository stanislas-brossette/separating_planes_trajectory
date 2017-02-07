#include <cstdlib>
#include <ctime>
#include <iostream>

#include <manifolds/CartesianProduct.h>

#include <pgsolver/utils/finiteDiff.h>
#include <pgsolver/utils/loadFromYaml.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverTrustRegionFilter.h>

#include <feet-trajectory/BoxTrajProblem.hh>
#include <feet-trajectory/utils/Printer.hh>

using namespace feettrajectory;

int main(void)
{
  std::cout << "StraightTraj Test: In this test we expect to generate a sequence of successive, non colliding, feet positions. It should end up describing a straight line" << std::endl;
  std::srand(static_cast<uint>(std::time(0)));
  std::string testDir = TESTS_DATA_DIR;
  std::string ymlPath = testDir + "/straightTraj.yml";
  ProblemConfig config(ymlPath);

  int nBoxes = config["nBoxes"];
  int nObstacles = config["nObstacles"];
  std::cout << "nBoxes: " << nBoxes << std::endl;

  mnf::CartesianProduct* M = BoxTrajProblemOnManifold::buildManifold(nBoxes, nObstacles);
  M->display();

  BoxTrajProblemOnManifold myProb(*M, ymlPath);

  Eigen::VectorXd v0(M->representationDim());
  v0 = myProb.findInitPoint();

  std::cout << "v0: " << v0 << std::endl;

  pgs::utils::finiteDiffCheck(myProb);

  myProb.M().forceOnM(v0, v0);
  mnf::Point x0 = myProb.M().createPoint(v0);
  std::cout << "x0 = \n" << x0 << std::endl;

  // Solver mySolver;
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  pgs::Results res = mySolver.solve();

  std::cout << "xSol = \n" << res.x_star << std::endl;

  print("log.yml", myProb, res.x_star);

  return 0;
}
