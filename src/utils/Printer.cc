#include <feet-trajectory/utils/Printer.hh>
#include <Eigen/Core>

namespace feettrajectory
{
void print(const std::string& fileName, const BoxTrajProblemOnManifold& pb,
           const mnf::Point& xStar)
{
  Eigen::IOFormat logFmt(3, 0, ", ", "\n", "[", "]");
  std::ofstream logFile;
  logFile.open(fileName);

  // Boxes
  auto bSize = pb.boxSize().transpose().format(logFmt);
  logFile << "InitialBox: {size: " << bSize
          << ", position: " << pb.initPos().transpose().format(logFmt) << "}\n";
  logFile << "MobileBoxes:\n";
  for (size_t i = 0; i < pb.nBoxes(); i++)
  {
    logFile << "- {size: " << bSize
            << ", position: " << xStar(0)(i).value().transpose().format(logFmt)
            << "}\n";
  }
  logFile << "FinalBox: {size: " << bSize
          << ", position: " << pb.finalPos().transpose().format(logFmt)
          << "}\n";

  // Obstacles
  logFile << "Obstacles:\n";
  for (size_t i = 0; i < pb.nObstacles(); i++)
  {
    logFile << "- {size: "
            << pb.obstacles()[i].size().transpose().format(logFmt)
            << ", position: "
            << pb.obstacles()[i].center().transpose().format(logFmt) << "}"
            << std::endl;
  }

  // Separating planes
  logFile << "SeparatingPlanes:\n";
  for (size_t i = 0; i < pb.nPlans(); i++)
    logFile << "- { d: " << xStar(1)(i)(0).value()
            << ", normal: " << xStar(1)(i)(1).value().transpose().format(logFmt)
            << ", boxAbove: " << pb.plans()[i].boxAbove()
            << ", obstacleBelow: " << pb.plans()[i].boxBelow() << "}"
            << std::endl;
}
} /* feettrajectory */
