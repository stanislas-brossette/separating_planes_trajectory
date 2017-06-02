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

int main(int argc, char* argv[])
{
  std::string ymlPath;
  if (argc > 1)
  {
    ymlPath = std::string(CONFIGS_DATA_DIR) + "/" + argv[1] + ".yml";
  }
  else
  {
    std::cout << "Loading default file \"singleObstacle.yml\"" << std::endl;
    ymlPath = std::string(CONFIGS_DATA_DIR) + "/singleObstacle.yml";
  }
  ProblemConfig config(ymlPath);

  int nBoxes = config["nBoxes"];
  int nObstacles = 0;
  if (config.has("obstacles"))
    nObstacles = static_cast<int>(config["obstacles"].asVecBox().size());

  mnf::CartesianProduct* M =
      BoxesHullTrajProblem::buildManifold(nBoxes, nObstacles);
  M->display();

  BoxesHullTrajProblem myProb(*M, ymlPath);

  // Initialization
  Eigen::VectorXd initVec(myProb.dimVar());
  if (myProb.config().has("x0"))
  {
    initVec << myProb.config()["x0"].asVectorXd();
    std::cout << "myProb.config()[x0].asVectorXd(): "
              << myProb.config()["x0"].asVectorXd().transpose() << std::endl;
  }
  else
    initVec = myProb.findInitPoint();

  if (myProb.config().has("initialGuessRandomFactor"))
  {
    initVec = initVec +
              myProb.config()["initialGuessRandomFactor"] *
                  Eigen::VectorXd::Random(myProb.dimVar());
  }

  myProb.M().forceOnM(initVec, initVec);

  pgs::utils::finiteDiffCheck(myProb);

  mnf::Point x0 = myProb.M().createPoint(initVec);

  // Solver mySolver;
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  pgs::Results res = mySolver.solve();

  if (config.has("plotResult") && config["plotResult"])
  {
    if (config.has("logLevel") && config["logLevel"].compare("NO_LOG") != 0)
    {
      std::cout << "Plotting all iterations" << std::endl;
      printAllIterations(config["logName"], myProb, res.x_star,
                         mySolver.logger().folder());
      std::string command = "./animSteps.py logs/";
      if (argc > 1)
      {
        command += argv[1];
      }
      else
      {
        command += "singleObstacle";
      }
      command += ".log";
      system(command.c_str());
    }
  }
  return 0;
}
