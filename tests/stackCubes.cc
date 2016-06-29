#include <cstdlib>
#include <ctime>

#include <pgsolver/utils/finiteDiff.h>
#include <pgsolver/utils/loadFromYaml.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverTrustRegionFilter.h>
#include <pgsolver/utils/TimeLogger.h>

#include <manifolds/CartesianProduct.h>
#include <manifolds/Point.h>

#include <cube-stacks/CubeStackProblemOnManifold.hh>
#include <cube-stacks/utils/ProblemConfig.hh>

using namespace cubestacks;

int main(void)
{
  std::srand(static_cast<uint>(std::time(0)));
  std::string testDir = TESTS_DATA_DIR;
  std::string ymlPath = testDir + "/stackCubesTmp.yml";
  ProblemConfig config(ymlPath);
  int nCubes = config["nCubes"];

  mnf::CartesianProduct* M = CubeStackProblemOnManifold::buildManifold(nCubes);

  M->display();

  CubeStackProblemOnManifold myProb(*M, ymlPath);
  Eigen::VectorXd v0(M->representationDim());
  v0 = myProb.findInitPoint();
  //M->createRandomPoint(v0);

  //for (size_t i = 0; i < myProb.nCubes_; i++)
  //{
    //Eigen::VectorXd vCube(7);
    //vCube.head(3) =
        //config["init.cube" + std::to_string(i) + ".trans"].asVector3d();
    //vCube.segment(3, 4) =
        //config["init.cube" + std::to_string(i) + ".quat"].asVector4d();
    //v0.segment(static_cast<long>(7 * i), 7) << vCube;
  //}
  //for (size_t i = 0; i < myProb.nPlans_; i++)
  //{
    //Eigen::VectorXd vPlan(4);
    //vPlan.head(1) << config["init.plan" + std::to_string(i) + ".d"].asDouble();
    //vPlan.segment(1, 3) =
        //config["init.plan" + std::to_string(i) + ".normal"].asVector3d();
    //v0.segment(static_cast<long>(7 * myProb.nCubes_ + 4 * i), 4) << vPlan;
  //}

  std::cout << "v0: " << v0 << std::endl;

  pgs::utils::finiteDiffCheck(myProb);

  myProb.M().forceOnM(v0, v0);
  mnf::Point x0 = myProb.M().createPoint(v0);
  std::cout << "x0 = \n" << x0 << std::endl;

  // Solver mySolver;
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  myProb.fileForMatlab(
      "/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      "stack" + std::to_string(nCubes) + "CubesInit.m",
      x0);
  pgs::Results res = mySolver.solve();
  //auto timeLog = mySolver.getTimeLogger();
  //std::cout << *timeLog << std::endl;

  std::cout << "xSol = \n" << res.x_star << std::endl;
  myProb.fileForMatlab(
      "/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      "stack" + std::to_string(nCubes) + "Cubes.m",
      res.x_star);
  return 0;
}
