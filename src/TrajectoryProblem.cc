#include <iostream>
#include <fstream>
#include <limits>

#include <Eigen/Geometry>

#include <feet-trajectory/TrajectoryProblem.hh>

namespace feettrajectory
{
TrajectoryProblem::TrajectoryProblem(const std::string& configPath)
    : config_(configPath),
      initPos_(config_["initPos"].asVector3d()),
      finalPos_(config_["finalPos"].asVector3d()),
      boxSize_(config_["BoxSize"].asVector3d()),
      nBoxes_(config_["nBoxes"].asSize_t()),
      initBox_(-1, boxSize_),
      finalBox_(static_cast<int>(nBoxes_), boxSize_),
      costFct_(nBoxes_, initPos_, finalPos_),
      initBoxAbovePlanFct_(initBox_),
      finalBoxAbovePlanFct_(finalBox_)
{
  if (config_.has("obstacles")) obstacles_ = config_["obstacles"].asVecBox();

  if (config_.has("fixedPlanes"))
    fixedPlanes_ = config_["fixedPlanes"].asVecFixedPlan();

  nObstacles_ = obstacles_.size();
  nFixedPlanes_ = fixedPlanes_.size();
  nPlans_ = (nBoxes_ + 1) * nObstacles_;

  nMobilePlanCstr_ = nObstacles_ * (nBoxes_ + 1);
  nFixedPlanCstr_ = nFixedPlanes_ * nBoxes_;
  numberOfCstr_ = 3 * nMobilePlanCstr_ + nFixedPlanCstr_;

  for (int i = 0; i < static_cast<int>(nBoxes_); ++i)
  {
    Box b(i, boxSize_);
    boxes_.push_back(b);
    boxAbovePlanFcts_.push_back(BoxAbovePlan(b));
  }

  for (auto p : fixedPlanes_)
  {
    for (auto b : boxes_)
    {
      boxAboveFixedPlanFcts_.push_back(BoxAboveFixedPlan(b, p));
    }
  }

  for (auto o : obstacles_)
  {
    obstacleAbovePlanFcts_.push_back(BoxAbovePlan(o));
  }

  for (int i = -1; i < static_cast<int>(nBoxes_); ++i)
  {
    for (int j = 0; j < static_cast<int>(nObstacles_); ++j)
    {
      plans_.push_back(PlanForHull(j, i, i + 1));
    }
  }

  for (auto c : boxAboveFixedPlanFcts_)
  {
    Eigen::IOFormat CleanFmt(2, 0, ", ", "", "[", "]");
    std::stringstream ss;
    ss << "Box " << std::to_string(c.box().index())
       << " above fixed plan: {normal:"
       << c.normal().transpose().format(CleanFmt) << ", d: " << c.d() << "}";
    std::cout << ss.str() << std::endl;
    cstrNames_.push_back(ss.str());
  }

  for (size_t i = 0; i < static_cast<size_t>(nPlans_); i++)
  {
    std::string str("Plan " + std::to_string(i) + " between HullBoxes " +
                    std::to_string(plans_[i].box0Above()) + " and " +
                    std::to_string(plans_[i].box1Above()) + " and obstacle " +
                    std::to_string(plans_[i].boxBelow()));
    std::cout << str << std::endl;
    cstrNames_.push_back(str);
  }
}

TrajectoryProblem::~TrajectoryProblem() {}

std::string TrajectoryProblem::getCstrName(const size_t) const
{
  std::string str("");
  return str;
}

std::ostream& TrajectoryProblem::print(std::ostream& o) const
{
  return o << "TrajectoryProblem: ";
}

std::ostream& operator<<(std::ostream& o, const TrajectoryProblem& f)
{
  return f.print(o);
}

Eigen::VectorXd TrajectoryProblem::getBoxPositionsFromX(
    const Eigen::VectorXd& x) const
{
  Eigen::VectorXd boxPos(3 * nBoxes_);
  boxPos << x.head(3 * nBoxes_);
  return boxPos;
}

Eigen::VectorXd TrajectoryProblem::getPlansNormalsFromX(
    const Eigen::VectorXd& x) const
{
  Eigen::VectorXd plansN(3 * nPlans_);
  for (long i = 0; i < nPlans_; i++)
  {
    plansN.segment<3>(3 * i) << x.segment<3>(3 * nBoxes_ + 4 * i + 1);
  }
  return plansN;
}

Eigen::VectorXd TrajectoryProblem::getPlansDistancesFromX(
    const Eigen::VectorXd& x) const
{
  Eigen::VectorXd plansD(nPlans_);
  for (long i = 0; i < nPlans_; i++)
  {
    plansD[i] = x[3 * nBoxes_ + 4 * i];
  }
  return plansD;
}

void TrajectoryProblem::normalizeNormals(Eigen::Ref<Eigen::VectorXd> x) const
{
  for (int i = 0; i < nPlans_; i++)
  {
    x.segment<3>(3 * nBoxes_ + 4 * i + 1).normalize();
  }
}
} /* feettrajectory */
