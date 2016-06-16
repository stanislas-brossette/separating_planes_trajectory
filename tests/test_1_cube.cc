//#include <cube-stacks/CubeStackProblem.hh>
//#include <cube-stacks/CubeStackProblem1Cube.hh>

#include <manifolds/CartesianProduct.h>
#include <cube-stacks/CubeStackProblemOnManifold.hh>
#include <cube-stacks/utils/ProblemConfig.hh>

using namespace cubestacks;

int main(void)
{
  std::string testDir = TESTS_DATA_DIR;
  ProblemConfig config(testDir + "/test.yml");

  mnf::CartesianProduct* M = CubeStackProblemOnManifold::buildManifold(1);

  CubeStackProblemOnManifold test( *M, testDir + "/test.yml");
  //CubeStackProblem1Cube test1;
  Eigen::VectorXd x0(7);
  x0 << 1, 0, 10, 0, 10, 10, 10;
  //test.init(x0);
  //test.solve("pgsolver");
  return 0;
}
