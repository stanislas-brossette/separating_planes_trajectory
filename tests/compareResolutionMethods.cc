#include <iostream>

#include <manifolds/CartesianProduct.h>

#include <pgsolver/utils/finiteDiff.h>
#include <pgsolver/utils/loadFromYaml.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverTrustRegionFilter.h>

#include <feet-trajectory/BoxesHullTrajProblem.hh>
#include <feet-trajectory/TrajectoryProblem.hh>
#include <feet-trajectory/utils/Printer.hh>
#include <feet-trajectory/utils/ProblemConfig.hh>
#include <feet-trajectory/utils/AlternateQPSolver.hh>

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
    std::cout << "Loading default file \"trajProblem.yml\"" << std::endl;
    ymlPath = std::string(CONFIGS_DATA_DIR) + "/singleObstacle.yml";
  }

  ProblemConfig config(ymlPath);

  Eigen::VectorXd initVec;
  {
    TrajectoryProblem myProb(ymlPath);
    initVec.resize(myProb.dimVar());
    if (config.has("x0"))
    {
      initVec << config["x0"].asVectorXd();
      std::cout << "config[x0].asVectorXd(): "
                << config["x0"].asVectorXd().transpose() << std::endl;
    }
    else
      initVec = myProb.findInitPoint();

    if (config.has("initialGuessRandomFactor"))
    {
      initVec = initVec +
                config["initialGuessRandomFactor"] *
                    Eigen::VectorXd::Random(myProb.dimVar());
    }

    std::cout << "myProb: " << myProb << std::endl;
    AlternateQPSolver altQP(myProb, myProb.maxIter());

    myProb.normalizeNormals(initVec);

    altQP.init(initVec);

    altQP.solve();
    initVec = altQP.res();
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
  }
  
  {
    int nBoxes = config["nBoxes"];
    int nObstacles = 0;
    if (config.has("obstacles"))
      nObstacles = static_cast<int>(config["obstacles"].asVecBox().size());

    mnf::CartesianProduct* M =
        BoxesHullTrajProblem::buildManifold(nBoxes, nObstacles);
    M->display();

    BoxesHullTrajProblem myProb(*M, ymlPath);

    Eigen::VectorXd v0(M->representationDim());
    //v0 = myProb.findInitPoint();
    v0 = initVec;
    myProb.M().forceOnM(v0, v0);

    std::cout << "v0: " << v0.format(feettrajectory::fmt::custom) << std::endl;

    // pgs::utils::finiteDiffCheck(myProb);

    myProb.M().forceOnM(v0, v0);
    mnf::Point x0 = myProb.M().createPoint(v0);

    // Solver mySolver;
    pgs::SolverTrustRegionFilter mySolver;

    pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

    mySolver.init(myProb, x0);
    pgs::Results res = mySolver.solve();

    std::cout << "xSol = \n" << res.x_star << std::endl;

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

  return 0;
}
