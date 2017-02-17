#include <iostream>
#include <feet-trajectory/utils/Printer.hh>
#include <feet-trajectory/BoxTrajProblem.hh>
#include <feet-trajectory/BoxesHullTrajProblem.hh>
#include <Eigen/Core>

namespace feettrajectory
{
void print(const std::string& fileName, const BoxTrajProblem& pb,
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

void print(const std::string& fileName, const BoxesHullTrajProblem& pb,
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
            << ", box0Above: " << pb.plans()[i].box0Above()
            << ", box1Above: " << pb.plans()[i].box1Above()
            << ", obstacleBelow: " << pb.plans()[i].boxBelow() << "}"
            << std::endl;

  // Fixed Planes
  logFile << "FixedPlanes:\n";
  for (size_t i = 0; i < pb.nFixedPlanes(); i++)
  {
    logFile << "- { d: " << pb.fixedPlanes()[i].d() << ", normal: "
            << pb.fixedPlanes()[i].normal().transpose().format(logFmt) << "}"
            << std::endl;
  }
}

void print(const std::string& fileName, const Eigen::Vector3d& bSize,
           const Eigen::Vector3d& oSize, const Eigen::Vector3d& oPos,
           const Eigen::Vector3d& t, const double& d, const Eigen::Vector3d& n)
{
  Eigen::IOFormat logFmt(3, 0, ", ", "\n", "[", "]");
  std::ofstream logFile;
  logFile.open(fileName);

  // Boxes
  logFile << "MobileBoxes:\n";
  logFile << "- {size: " << bSize.transpose().format(logFmt)
          << ", position: " << t.transpose().format(logFmt) << "}\n";

  // Obstacles
  logFile << "Obstacles:\n";
  logFile << "- {size: " << oSize.transpose().format(logFmt)
          << ", position: " << oPos.transpose().format(logFmt) << "}"
          << std::endl;

  // Separating planes
  logFile << "SeparatingPlanes:\n";
  logFile << "- { d: " << d << ", normal: " << n.transpose().format(logFmt)
          << ", boxAbove: " << 0 << ", obstacleBelow: " << 0 << "}"
          << std::endl;
}

void printAllIterations(const std::string& fileName,
                        const BoxesHullTrajProblem& pb, const mnf::Point& xStar,
                        const std::string& folder)
{
  Eigen::IOFormat logFmt(3, 0, ", ", "\n", "[", "]");
  std::ofstream logFile;
  logFile.open(fileName);

  // Boxes
  auto bSize = pb.boxSize().transpose().format(logFmt);
  logFile << "InitialBox: {size: " << bSize
          << ", position: " << pb.initPos().transpose().format(logFmt) << "}\n";
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
  // Fixed Planes
  logFile << "FixedPlanes:\n";
  for (size_t i = 0; i < pb.nFixedPlanes(); i++)
  {
    logFile << "- { d: " << pb.fixedPlanes()[i].d() << ", normal: "
            << pb.fixedPlanes()[i].normal().transpose().format(logFmt) << "}"
            << std::endl;
  }

  std::string xFile(folder + "xLog.m");
  std::vector<Eigen::VectorXd> all_x = parseX(xFile, xStar);
  logFile << "nIter: " << all_x.size() << std::endl;
  for (size_t iter = 0; iter < all_x.size(); iter++)
  {
    mnf::Point xIter = xStar.getManifold().createPoint(all_x[iter]);
    // Separating Planes
    logFile << "SeparatingPlanes" << iter << ":\n";
    for (size_t i = 0; i < pb.nPlans(); i++)
      logFile << "- { d: " << xIter(1)(i)(0).value() << ", normal: "
              << xIter(1)(i)(1).value().transpose().format(logFmt)
              << ", box0Above: " << pb.plans()[i].box0Above()
              << ", box1Above: " << pb.plans()[i].box1Above()
              << ", obstacleBelow: " << pb.plans()[i].boxBelow() << "}"
              << std::endl;
    // Mobile Boxes
    logFile << "MobileBoxes" << iter << ":\n";
    for (size_t i = 0; i < pb.nBoxes(); i++)
    {
      logFile << "- {size: " << bSize << ", position: "
              << xIter(0)(i).value().transpose().format(logFmt) << "}\n";
    }
  }
}

std::vector<Eigen::VectorXd> parseX(const std::string& fileName,
                                    const mnf::Point& x)
{
  std::vector<Eigen::VectorXd> res;
  std::ifstream file(fileName, std::ios::in);
  std::string s;
  while (std::getline(file, s))
  {
    if (s.compare(0, 2, "x_") == 0)
    {
      Eigen::VectorXd newX(x.getRepresentationDimM());
      size_t posNumber = s.find("= [[");
      size_t posNumberEnd = s.find("]");
      std::string subStr =
          s.substr(posNumber + 5, posNumberEnd - posNumber - 5);
      std::string::iterator end_pos =
          std::remove(subStr.begin(), subStr.end(), ' ');
      subStr.erase(end_pos, subStr.end());
      newX[0] = std::stod(subStr);
      for (long i = 1; i < newX.size(); i++)
      {
        std::getline(file, s);
        posNumber = s.find("[");
        posNumberEnd = s.find("]");
        subStr = s.substr(posNumber + 1, posNumberEnd - posNumber - 1);
        newX[i] = std::stod(subStr);
      }
      res.push_back(newX);
    }
  }
  return res;
}

} /* feettrajectory */
