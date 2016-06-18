
#include <pgsolver/utils/finiteDiff.h>
#include <pgsolver/utils/loadFromYaml.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverTrustRegionFilter.h>

#include <manifolds/CartesianProduct.h>
#include <manifolds/Point.h>

#include <cube-stacks/CubeStackProblemOnManifold.hh>
#include <cube-stacks/utils/ProblemConfig.hh>

using namespace cubestacks;

int main(void)
{
  std::string testDir = TESTS_DATA_DIR;
  std::string ymlPath = testDir + "/test.yml";
  ProblemConfig config(ymlPath);
  int nCubes = config["nCubes"];

  mnf::CartesianProduct* M = CubeStackProblemOnManifold::buildManifold(nCubes);

  M->display();

  CubeStackProblemOnManifold myProb( *M, testDir + "/test.yml");
  //CubeStackProblem1Cube test1;
  Eigen::VectorXd v0(M->representationDim());
  v0 = myProb.findInitPoint();

  //for (int i = 0; i < nCubes; i++)
  //{
    //Eigen::VectorXd vCube(7);
    //vCube.head(3) = config["init.cube"+std::to_string(i)+".trans"].asVector3d();
    //vCube.segment(3,4) = config["init.cube"+std::to_string(i)+".quat"].asVector4d();
    //v0.segment(7*i,7) << vCube;
  //}


  std::cout << "v0: " << v0 << std::endl;

  pgs::utils::finiteDiffCheck(myProb);

  myProb.M().forceOnM(v0, v0);
  mnf::Point x0 = myProb.M().createPoint(v0);
  std::cout << "x0 = \n" << x0 << std::endl;

  //Solver mySolver;
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  pgs::Results res = mySolver.solve();
  //test.init(x0);
  //test.solve("pgsolver");
  return 0;
}
