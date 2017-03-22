#include <iostream>

#include <feet-trajectory/TrajectoryProblem.hh>
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
  AlternateQPSolver altQP(myProb);

  Eigen::VectorXd initVec(myProb.dimVar());
  std::cout << "myProb.config()[x0].asVectorXd(): "
            << myProb.config()["x0"].asVectorXd() << std::endl;
  if (myProb.config().has("x0"))
    initVec << myProb.config()["x0"].asVectorXd();
  else
    initVec.setRandom();

  myProb.normalizeNormals(initVec);

  altQP.init(initVec);

  std::cout << "altQP.QPWithNFixed(): " << altQP.qpNfixed() << std::endl;

  altQP.solve();

  return 0;
}
