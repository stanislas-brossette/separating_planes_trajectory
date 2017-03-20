#include <iostream>

#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/utils/AlternateQPSolver.hh>

using namespace feettrajectory;

int main()
{
  // assert(argc > 1 && "Please provide a file to load");
  // std::cout << "Generate Trajectory: Loads a config and generates a
  // trajectory "
  //"between initPos and finalPos that avoids the obstacles and "
  //"some fixed planes" << std::endl;

  // std::string ymlPath = std::string(CONFIGS_DATA_DIR) + "/" + argv[1] +
  // ".yml";
  std::string ymlPath = std::string(CONFIGS_DATA_DIR) + "/trajProblem.yml";

  TrajectoryProblem myProb(ymlPath);
  std::cout << "myProb: " << myProb << std::endl;
  AlternateQPSolver altQP(myProb);

  return 0;
}
