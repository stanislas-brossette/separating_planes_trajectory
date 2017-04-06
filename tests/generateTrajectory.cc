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
  assert(argc > 1 && "Please provide a file to load");
  std::cout << "Generate Trajectory: Loads a config and generates a trajectory "
               "between initPos and finalPos that avoids the obstacles and "
               "some fixed planes" << std::endl;

  std::srand(static_cast<uint>(std::time(0)));
  std::string ymlPath = std::string(CONFIGS_DATA_DIR) + "/" + argv[1] + ".yml";
  ProblemConfig config(ymlPath);

  int nBoxes = config["nBoxes"];
  int nObstacles = 0;
  if (config.has("obstacles"))
    nObstacles = static_cast<int>(config["obstacles"].asVecBox().size());

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

  //pgs::utils::finiteDiffCheck(myProb);

  myProb.M().forceOnM(v0, v0);
  mnf::Point x0 = myProb.M().createPoint(v0);

  // Solver mySolver;
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  pgs::Results res = mySolver.solve();

  std::cout << "xSol = \n" << res.x_star << std::endl;

  if (config.has("plotResult") && config["plotResult"])
  {
    if (config.has("logLevel") && config["logLevel"].compare("NO_LOG") != 0)
    {
      std::cout << "Plotting all iterations" << std::endl;
      printAllIterations(config["logName"], myProb, res.x_star,
                         mySolver.logger().folder());
      std::string command = "./animSteps.py logs/";
      command += argv[1];
      command += ".log";
      system(command.c_str());
      return 0;
    }
    else
    {
      print(config["logName"], myProb, res.x_star);
      // Call the Python script passing a filename argument.
      std::string command = "./viewSteps.py logs/";
      command += argv[1];
      command += ".log";
      system(command.c_str());
      return 0;
    }
  }

  return 0;
}
