#include <iostream>

#include <feet-trajectory/TrajectoryProblem.hh>
//#include <feet-trajectory/utils/Printer.hh>
#include <feet-trajectory/utils/ProblemConfig.hh>
#include <feet-trajectory/utils/AlternateQPSolverJerk.hh>
#include <feet-trajectory/utils/Jerk3DIntegrator.hh>

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

  Eigen::VectorXd initVecPos = initVec.head(myProb.dimBoxes());
  Jerk3DIntegrator integ(0.3, myProb.nBoxes());
  Eigen::VectorXd state0(9);
  state0.setZero();
  state0.segment(0, 3) = myProb.initPos();
  integ.jerkFromPos(initVec.head(myProb.dimBoxes()),
                    initVec.head(myProb.dimBoxes()), state0, 1e-3);

  std::cout << "myProb: " << myProb << std::endl;

  // Normalize the initial guess (to have correct normals)
  myProb.normalizeNormals(initVec);

  // Instanciate the solver
  AlternateQPSolverJerk altQP(myProb, myProb.maxIter(), integ, state0);

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
