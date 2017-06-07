#include <iostream>

#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/utils/Printer.hh>
#include <feet-trajectory/utils/ProblemConfig.hh>
#include <feet-trajectory/utils/AlternateQPSolver.hh>

using namespace feettrajectory;

int main(int argc, char* argv[])
{
  // Setup the problem
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
  TrajectoryProblem myProb(ymlPath);

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

  std::cout << "myProb: " << myProb << std::endl;

  // Instanciate the solver
  AlternateQPSolver altQP(myProb, myProb.maxIter());

  // Normalize the initial guess (to have correct normals)
  myProb.normalizeNormals(initVec);

  // Solve the problem
  altQP.init(initVec);
  altQP.solve();

  // Log the results
  altQP.logAllX("logs/altQP/");

  // Display the results
  printAllIterations(myProb.config()["logName"], myProb, altQP.res(),
                     "logs/altQP/");
  std::string command = "./animSteps.py logs/";
  if (argc > 1)
  {
    command += argv[1];
  }
  else
  {
    command += "singleObstacle";
  }
  command += ".log ";
  if (myProb.config().has("plotPlanes"))
    command += myProb.config()["plotPlanes"];
  system(command.c_str());

  return 0;
}
