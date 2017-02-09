#include <iostream>
#include <fstream>
#include <limits>

#include <Eigen/Geometry>

//#include <pgsolver/utils/usingManifold.h>

#include <manifolds/RealSpace.h>
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>

#include <feet-trajectory/BoxesHullTrajProblem.hh>
#include <feet-trajectory/functions/BoxAbovePlan.hh>

using mnf::RealSpace;
using mnf::SO3;
using mnf::S2;
using mnf::CartesianProduct;
using mnf::ExpMapMatrix;
using mnf::ExpMapQuaternion;
using mnf::Manifold;

namespace feettrajectory
{
BoxesHullTrajProblem::BoxesHullTrajProblem(const Manifold& M,
                                           const std::string& configPath)
    : Problem(M),
      config_(configPath),
      nBoxes_(config_["nBoxes"].asInt()),
      nObstacles_(config_["nObstacles"].asInt()),
      boxSize_(config_["BoxSize"].asVector3d()),
      initPos_(config_["initPos"].asVector3d()),
      finalPos_(config_["finalPos"].asVector3d()),
      initBox_(-1, boxSize_),
      finalBox_(nBoxes_, boxSize_),
      initBoxAbovePlanFct_(initBox_),
      finalBoxAbovePlanFct_(finalBox_)
{
  nPlans_ = (nBoxes_ + 1) * nObstacles_;

  for (int i = 0; i < static_cast<int>(nBoxes_); ++i)
  {
    Box b(i, boxSize_);
    boxes_.push_back(b);
    boxAbovePlanFcts_.push_back(BoxAbovePlan(b));
  }

  for (int i = 0; i < static_cast<int>(nObstacles_); ++i)
  {
    Box o(i, config_["ObstacleSizes"].asVector3d(),
          config_["ObstaclePositions"].asVector3d());
    obstacles_.push_back(o);
    obstacleAbovePlanFcts_.push_back(BoxAbovePlan(o));
  }

  for (int i = -1; i < static_cast<int>(nBoxes_); ++i)
  {
    for (int j = 0; j < nObstacles_; ++j)
    {
      plans_.push_back(PlanForHull(j, i, i + 1));
    }
  }

  for (int i = 0; i < nPlans_; i++)
  {
    std::string str("Plan" + std::to_string(i) + "BetweenHullBoxes" +
                    std::to_string(plans_[i].box0Above()) + "And" +
                    std::to_string(plans_[i].box1Above()) + "AndObstacle" +
                    std::to_string(plans_[i].boxBelow()));
    std::cout << str << std::endl;
    cstrNames_.push_back(str);
  }

  std::stringstream sstm;
  sstm << "FeetTraj" << nBoxes_ << "Boxes";
  name() = sstm.str();

  outRepObjDiff_.resize(1, M.representationDim());
  outRepObjDiff_.setZero();
  outRep_.resize(static_cast<Index>(24 * nPlans_), M.representationDim());
  outRep_.setZero();
}

BoxesHullTrajProblem::~BoxesHullTrajProblem() {}

Eigen::VectorXd BoxesHullTrajProblem::findInitPoint()
{
  mnf::Point xRand = M().createRandomPoint();
  // Eigen::Vector3d interBoxes;
  // for (size_t i = 0; i < nBoxes_; ++i)
  //{
  // xRand(0)(i)[0][2] += 2 * static_cast<double>(i);
  //}
  // for (size_t i = 0; i < nPlans_; ++i)
  //{
  // size_t boxBelow = plans_[i].boxBelow();
  // size_t boxAbove = plans_[i].boxAbove();
  // Eigen::Vector3d boxBelowPos = xRand(0)(boxBelow)[0];
  // Eigen::Vector3d boxAbovePos = xRand(0)(boxAbove)[0];
  // interBoxes = boxAbovePos - boxBelowPos;
  // interBoxes = interBoxes / interBoxes.norm();
  // xRand(1)(i)[1] = interBoxes;
  // xRand(1)(i)[0][0] = interBoxes.dot((boxAbovePos + boxBelowPos) / 2);
  //}
  return xRand.value();
}

CartesianProduct* BoxesHullTrajProblem::buildManifold(const Index& nBoxes,
                                                      const Index& nObstacles)
{
  assert(nBoxes > 0 && "Number of Boxes must be positive and non null");
  assert(nObstacles >= 0 && "Number of Obstacles must be positive");
  RealSpace* R1_ = new RealSpace(1);
  S2* S2_ = new S2();
  RealSpace* MBox = new RealSpace(3);
  CartesianProduct* MPlane = new CartesianProduct(*R1_, *S2_);
  RealSpace* R0 = new RealSpace(0);
  CartesianProduct* MBoxes;
  if (nBoxes == 1)
    MBoxes = new CartesianProduct(*MBox, *R0);
  else
  {
    MBoxes = new CartesianProduct(*MBox, *MBox);
    for (size_t i = 2; i < static_cast<size_t>(nBoxes); ++i)
      MBoxes->multiply(*MBox);
  }

  Index nPlanes = (nBoxes + 1) * nObstacles;
  CartesianProduct* MPlanes;
  if (nPlanes == 0)
    MPlanes = new CartesianProduct(*R0, *R0);
  else if (nPlanes == 1)
    MPlanes = new CartesianProduct(*MPlane, *R0);
  else
  {
    MPlanes = new CartesianProduct(*MPlane, *MPlane);
    for (size_t i = 2; i < static_cast<size_t>(nPlanes); ++i)
      MPlanes->multiply(*MPlane);
  }
  CartesianProduct* MBoxesAndPlanes = new CartesianProduct(*MBoxes, *MPlanes);
  return MBoxesAndPlanes;
}

void BoxesHullTrajProblem::getTangentLB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = -infinity;
}

void BoxesHullTrajProblem::getTangentUB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = infinity;
}

void BoxesHullTrajProblem::evalObj(double& out) const
{
  out = 0;
  // distance between first mobile box and initial position
  Eigen::Vector3d pos = initPos_;
  Eigen::Vector3d posNext = phi_x_z()(0)(0)[0];
  Eigen::Vector3d dist = posNext - pos;
  out += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
  for (size_t i = 0; i < nBoxes_ - 1; ++i)
  {
    // distance between successive mobile boxes
    pos = phi_x_z()(0)(i)[0];
    posNext = phi_x_z()(0)(i + 1)[0];
    dist = posNext - pos;
    out += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
  }
  // distance between last mobile box and final position
  pos = phi_x_z()(0)(nBoxes_ - 1)[0];
  posNext = finalPos_;
  dist = posNext - pos;
  out += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
}
void BoxesHullTrajProblem::evalObjDiff(RefMat out) const
{
  assert(out.rows() == 1 && "wrong rows size");
  assert(out.cols() == M().tangentDim() && "wrong cols size");
  outRepObjDiff_.setZero();

  Index boxRepDim = M()(0)(static_cast<size_t>(0)).representationDim();
  Eigen::Vector3d pos = initPos_;
  Eigen::Vector3d posNext = phi_x_z()(0)(0)[0];
  Eigen::Vector3d dist = posNext - pos;
  outRepObjDiff_.block(0, (0) * boxRepDim, 1, boxRepDim) +=
      2 * dist.transpose();
  for (Index i = 0; i < static_cast<Index>(nBoxes_) - 1; ++i)
  {
    pos = phi_x_z()(0)(static_cast<size_t>(i))[0];
    posNext = phi_x_z()(0)(static_cast<size_t>(i) + 1)[0];
    dist = posNext - pos;
    outRepObjDiff_.block(0, i * boxRepDim, 1, boxRepDim) +=
        -2 * dist.transpose();
    outRepObjDiff_.block(0, (i + 1) * boxRepDim, 1, boxRepDim) +=
        2 * dist.transpose();
  }
  pos = phi_x_z()(0)(nBoxes_ - 1)[0];
  posNext = finalPos_;
  dist = posNext - pos;
  outRepObjDiff_.block(0, (static_cast<Index>(nBoxes_) - 1) * boxRepDim, 1,
                       boxRepDim) += -2 * dist.transpose();
  M().applyDiffRetractation(out, outRepObjDiff_, x().value());
}

void BoxesHullTrajProblem::evalLinCstr(RefVec, size_t) const {}

void BoxesHullTrajProblem::evalLinCstrDiff(RefMat, size_t) const {}

void BoxesHullTrajProblem::getLinCstrLB(RefVec, size_t) const {}

void BoxesHullTrajProblem::getLinCstrUB(RefVec, size_t) const {}

void BoxesHullTrajProblem::evalNonLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  out.setZero();

  auto iPlan = i;
  auto iBox0Above = plans_[iPlan].box0Above();
  auto iBox1Above = plans_[iPlan].box1Above();
  auto iBoxBelow = plans_[iPlan].boxBelow();
  Eigen::Vector4d nullQuat(0, 0, 0, 1);
  Eigen::Vector3d transBelow = obstacles_[iBoxBelow].center();
  Eigen::Vector3d normal = phi_x_z()(1)(iPlan)[1];
  double d = phi_x_z()(1)(iPlan)[0](0);

  Eigen::Vector3d trans0Above, trans1Above;

  const BoxAbovePlan* box0AbovePlanFct = &boxAbovePlanFcts_[iBox0Above];
  const BoxAbovePlan* box1AbovePlanFct = &boxAbovePlanFcts_[iBox1Above];

  if (iBox0Above != -1 && iBox1Above != nBoxes_)
  {
    trans0Above = phi_x_z()(0)(iBox0Above)[0];
    trans1Above = phi_x_z()(0)(iBox1Above)[0];
  }
  else if (iBox0Above == -1 && iBox1Above != nBoxes_)
  {
    trans0Above = initPos_;
    trans1Above = phi_x_z()(0)(iBox1Above)[0];
    box0AbovePlanFct = &initBoxAbovePlanFct_;
  }
  else if (iBox0Above != -1 && iBox1Above == nBoxes_)
  {
    trans0Above = phi_x_z()(0)(iBox0Above)[0];
    trans1Above = finalPos_;
    box1AbovePlanFct = &finalBoxAbovePlanFct_;
  }
  else
  {
    std::cerr << "Box0 is initial pos and Box1 is final, there is no "
                 "intermediary boxes..." << std::endl;
  }

  box0AbovePlanFct->compute(out.head(8), trans0Above, nullQuat, d, normal);
  box1AbovePlanFct->compute(out.segment(8, 8), trans1Above, nullQuat, d,
                            normal);
  obstacleAbovePlanFcts_[iBoxBelow].compute(out.tail(8), transBelow, nullQuat,
                                            -d, -normal);
}

void BoxesHullTrajProblem::evalNonLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == nonLinCstrDim(i) && "wrong row size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");

  auto iPlan = i;
  auto iBox0Above = plans_[iPlan].box0Above();
  auto iBox1Above = plans_[iPlan].box1Above();
  auto iBoxBelow = plans_[iPlan].boxBelow();

  Eigen::Vector4d nullQuat(0, 0, 0, 1);
  Eigen::Vector3d transBelow = obstacles_[iBoxBelow].center();
  Eigen::Vector3d normal = phi_x_z()(1)(iPlan)[1];
  double d = phi_x_z()(1)(iPlan)[0](0);
  auto rowBegin = 24 * iPlan;
  auto box0AboveBegin = 3 * iBox0Above;
  auto box1AboveBegin = 3 * iBox1Above;
  auto planBegin = 3 * nBoxes_ + 4 * iPlan;
  long rowBeginLong = static_cast<long>(rowBegin);
  long box0AboveBeginLong = static_cast<long>(box0AboveBegin);
  long box1AboveBeginLong = static_cast<long>(box1AboveBegin);
  long planBeginLong = static_cast<long>(planBegin);

  Eigen::Matrix<double, 8, 1> tmpWTF;

  Eigen::Vector3d trans0Above, trans1Above;
  const BoxAbovePlan* box0AbovePlanFct = &boxAbovePlanFcts_[iBox0Above];
  const BoxAbovePlan* box1AbovePlanFct = &boxAbovePlanFcts_[iBox1Above];
  if (iBox0Above != -1 && iBox1Above != nBoxes_)
  {
    trans0Above = phi_x_z()(0)(iBox0Above)[0];
    trans1Above = phi_x_z()(0)(iBox1Above)[0];
    box0AbovePlanFct->diffTrans(
        outRep_.block(rowBeginLong, box0AboveBeginLong, 8, 3), trans0Above,
        nullQuat, d, normal);
    box1AbovePlanFct->diffTrans(
        outRep_.block(rowBeginLong + 8, box1AboveBeginLong, 8, 3), trans1Above,
        nullQuat, d, normal);
  }
  else if (iBox0Above == -1 && iBox1Above != nBoxes_)
  {
    trans0Above = initPos_;
    trans1Above = phi_x_z()(0)(iBox1Above)[0];
    box0AbovePlanFct = &initBoxAbovePlanFct_;
    box1AbovePlanFct->diffTrans(
        outRep_.block(rowBeginLong + 8, box1AboveBeginLong, 8, 3), trans1Above,
        nullQuat, d, normal);
  }
  else if (iBox0Above != -1 && iBox1Above == nBoxes_)
  {
    trans0Above = phi_x_z()(0)(iBox0Above)[0];
    trans1Above = finalPos_;
    box1AbovePlanFct = &finalBoxAbovePlanFct_;
    box0AbovePlanFct->diffTrans(
        outRep_.block(rowBeginLong, box0AboveBeginLong, 8, 3), trans0Above,
        nullQuat, d, normal);
  }
  else
  {
    std::cerr << "Box0 is initial pos and Box1 is final, there is no "
                 "intermediary boxes..." << std::endl;
  }

  // Box0
  box0AbovePlanFct->diffD(tmpWTF, trans0Above, nullQuat, d, normal);
  outRep_.block(rowBeginLong, planBeginLong, 8, 1) = tmpWTF;
  box0AbovePlanFct->diffNormal(
      outRep_.block(rowBeginLong, planBeginLong + 1, 8, 3), trans0Above,
      nullQuat, d, normal);

  // Box1
  box1AbovePlanFct->diffD(tmpWTF, trans1Above, nullQuat, d, normal);
  outRep_.block(rowBeginLong + 8, planBeginLong, 8, 1) = tmpWTF;
  box1AbovePlanFct->diffNormal(
      outRep_.block(rowBeginLong + 8, planBeginLong + 1, 8, 3), trans1Above,
      nullQuat, d, normal);

  // Obstacle
  obstacleAbovePlanFcts_[iBoxBelow].diffD(tmpWTF, transBelow, nullQuat, -d,
                                          -normal);
  outRep_.block(rowBeginLong + 16, planBeginLong, 8, 1) = -tmpWTF;

  Eigen::Matrix<double, 8, 3> tmpDiffNormal;
  obstacleAbovePlanFcts_[iBoxBelow].diffNormal(tmpDiffNormal, transBelow,
                                               nullQuat, -d, -normal);
  outRep_.block(rowBeginLong + 16, planBeginLong + 1, 8, 3) = -tmpDiffNormal;

  M().applyDiffRetractation(out, outRep_.middleRows(rowBeginLong, 24),
                            x().value());

}

void BoxesHullTrajProblem::getNonLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  auto iPlan = i;
  auto iBox0Above = plans_[iPlan].box0Above();
  auto iBox1Above = plans_[iPlan].box1Above();
  auto iBoxBelow = plans_[iPlan].boxBelow();
  boxAbovePlanFcts_[iBox0Above].LB(out.head(8));
  boxAbovePlanFcts_[iBox1Above].LB(out.segment(8, 8));
  boxAbovePlanFcts_[iBoxBelow].LB(out.tail(8));
}
void BoxesHullTrajProblem::getNonLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  auto iPlan = i;
  auto iBox0Above = plans_[iPlan].box0Above();
  auto iBox1Above = plans_[iPlan].box1Above();
  auto iBoxBelow = plans_[iPlan].boxBelow();
  boxAbovePlanFcts_[iBox0Above].UB(out.head(8));
  boxAbovePlanFcts_[iBox1Above].UB(out.segment(8, 8));
  boxAbovePlanFcts_[iBoxBelow].UB(out.tail(8));
}

size_t BoxesHullTrajProblem::numberOfCstr() const { return nPlans_; }

Index BoxesHullTrajProblem::linCstrDim(size_t) const { return 0; }

Index BoxesHullTrajProblem::nonLinCstrDim(size_t) const { return 24; }

std::string BoxesHullTrajProblem::getCstrName(const size_t) const
{
  std::string str("");
  return str;
}

} /* feettrajectory */
