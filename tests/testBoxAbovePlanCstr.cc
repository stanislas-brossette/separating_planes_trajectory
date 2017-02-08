
// t: 0.304279 0.372069 0.720529
// q: 0 0 0 1
// d: 0.098299
// n: -0.586787 -0.381789  0.714085

#include <iostream>
#include <feet-trajectory/utils/ProblemConfig.hh>
#include <feet-trajectory/utils/Printer.hh>
#include <feet-trajectory/functions/BoxAbovePlan.hh>

using namespace feettrajectory;

int main(void)
{
  std::cout << "Test BoxAbovePlan constraint for various configurations"
            << std::endl;
  std::string testDir = TESTS_DATA_DIR;
  std::string ymlPath = testDir + "/testBoxAbovePlanCstr.yml";
  ProblemConfig config(ymlPath);

  Eigen::Vector3d boxPos(0., 0., 0.);
  Eigen::Vector3d boxSize(0.2, 0.1, 0.05);
  Eigen::Vector3d obsPos(0.5, 0.5, 0.5);
  Eigen::Vector3d obsSize(0.3, 1.2, 0.1);

  std::string testNumber = std::to_string(config["testNumber"].asInt());

  Eigen::Vector3d t = config["t" + testNumber].asVector3d();
  double planD = config["d" + testNumber].asDouble();
  Eigen::Vector3d planNormal = config["n" + testNumber].asVector3d();
  Eigen::Vector4d quat(0, 0, 0, 1);

  Box b(0, boxSize, boxPos);
  Box o(0, obsSize, obsPos);
  BoxAbovePlan cstrBox(b);
  BoxAbovePlan cstrObstacle(o);

  Eigen::VectorXd res(16);
  cstrBox.compute(res.head(8), t, quat, planD, planNormal);
  cstrObstacle.compute(res.tail(8), obsPos, quat, -planD, -planNormal);
  std::cout << "res: " << res << std::endl;

  print(config["logName"], boxSize, obsSize, obsPos, t, planD, planNormal);

  return 0;
}

