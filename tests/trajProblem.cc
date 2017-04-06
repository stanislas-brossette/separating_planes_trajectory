#include <iostream>

#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/utils/Printer.hh>
#include <feet-trajectory/utils/AlternateQPSolver.hh>

using namespace feettrajectory;

int main(int argc, char *argv[])
{
  std::cout << "Generate Trajectory: Loads a config and generates a trajectory"
               "between initPos and finalPos that avoids the obstacles and "
               "some fixed planes" << std::endl;

  std::string ymlPath;
  if (argc > 1)
  {
    ymlPath = std::string(CONFIGS_DATA_DIR) + "/" + argv[1] + ".yml";
  }
  else
  {
    std::cout << "Loading default file \"trajProblem.yml\"" << std::endl;
    ymlPath = std::string(CONFIGS_DATA_DIR) + "/singleObstacle.yml";
  }

  TrajectoryProblem myProb(ymlPath);
  std::cout << "myProb: " << myProb << std::endl;
  AlternateQPSolver altQP(myProb, myProb.maxIter());

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

  //initVec = Eigen::VectorXd::Random(myProb.dimVar());

  myProb.normalizeNormals(initVec);

  altQP.init(initVec);
  //std::cout << "qpPlanesFixed: " << altQP.qpPlanesFixed() << std::endl;
  //std::cout << "qpBoxesFixed: " << altQP.qpBoxesFixed() << std::endl;

  altQP.solve();
  altQP.logAllX("logs/altQP/");

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

  command += ".log";
  system(command.c_str());

  return 0;
}
