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
      maxIter_(config_["maxIter"].asSize_t()),
      initBox_(-1, boxSize_, initPos_, true),
      finalBox_(static_cast<int>(nBoxes_), boxSize_, finalPos_, true),
      initBoxAbovePlanFct_(initBox_),
      finalBoxAbovePlanFct_(finalBox_),
      costFct_(static_cast<long>(nBoxes_), initPos_, finalPos_)
{
  if (config_.has("obstacles"))
  {
    obstacles_ = config_["obstacles"].asVecBox();
    for (size_t i = 0; i < obstacles_.size(); i++)
    {
      obstacles_[i].setFixed(true);
    }
  }

  if (config_.has("fixedPlanes"))
    fixedPlanes_ = config_["fixedPlanes"].asVecFixedPlan();

  nObstacles_ = obstacles_.size();
  nFixedPlanes_ = fixedPlanes_.size();
  nPlans_ = (nBoxes_ + 1) * nObstacles_;

  nMobilePlanCstr_ = nObstacles_ * (nBoxes_ + 1);
  nFixedPlanCstr_ = nFixedPlanes_ * nBoxes_;
  numberOfCstr_ = 3 * nMobilePlanCstr_ + nFixedPlanCstr_;

  dimBox_ = 3;
  dimPlan_ = 4;
  dimNormal_ = 3;
  dimDistance_ = 1;
  dimBoxes_ = dimBox_ * nBoxes_;
  dimPlans_ = dimPlan_ * nPlans_;
  dimNormals_ = dimNormal_ * nPlans_;
  dimDistances_ = dimDistance_ * nPlans_;
  dimVar_ = dimBoxes_ + dimPlans_;

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

Eigen::VectorXd TrajectoryProblem::findInitPoint() const
{
  Eigen::VectorXd xInit(dimVar());
  xInit.setRandom();
  for (int i = 0; i < nBoxes(); i++)
  {
    xInit.segment<3>(3 * i) =
        initPos_ + (i + 1) * (finalPos_ - initPos_) / (nBoxes_ + 1);
  }
  return xInit;
}

Eigen::VectorXd TrajectoryProblem::getBoxPositionsFromX(
    const Eigen::VectorXd& x) const
{
  Eigen::VectorXd boxPos(3 * static_cast<long>(nBoxes_));
  boxPos << x.head(3 * static_cast<long>(nBoxes_));
  return boxPos;
}

Eigen::VectorXd TrajectoryProblem::getPlansNormalsFromX(
    const Eigen::VectorXd& x) const
{
  Eigen::VectorXd plansN(3 * static_cast<long>(nPlans_));
  for (long i = 0; i < nPlans_; i++)
  {
    plansN.segment<3>(3 * i)
        << x.segment<3>(3 * static_cast<long>(nBoxes_) + 4 * i + 1);
  }
  return plansN;
}

Eigen::VectorXd TrajectoryProblem::getPlansDistancesFromX(
    const Eigen::VectorXd& x) const
{
  Eigen::VectorXd plansD(nPlans_);
  for (long i = 0; i < nPlans_; i++)
  {
    plansD[i] = x[3 * static_cast<double>(nBoxes_) + 4 * i];
  }
  return plansD;
}

Eigen::Vector3d TrajectoryProblem::getBoxPositionFromX(
    size_t i, const Eigen::VectorXd& x) const
{
  return x.segment(3 * i, 3);
}

Eigen::Vector3d TrajectoryProblem::getPlanNormalFromX(
    size_t i, const Eigen::VectorXd& x) const
{
  return x.segment(dimBoxes_ + 4 * i + 1, 3);
}

double TrajectoryProblem::getPlanDistanceFromX(size_t i,
                                               const Eigen::VectorXd& x) const
{
  return x(dimBoxes_ + 4 * i);
}

void TrajectoryProblem::normalizeNormals(Eigen::Ref<Eigen::VectorXd> x) const
{
  for (int i = 0; i < nPlans_; i++)
  {
    double norm = x.segment<3>(dimBoxes_ + 4 * i + 1).norm();
    x.segment<3>(dimBoxes_ + 4 * i + 1).normalize();
    x(dimBoxes_ + 4 * i) /= norm;
  }
}
} /* feettrajectory */
