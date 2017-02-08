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

#include <feet-trajectory/BoxTrajProblem.hh>
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
BoxTrajProblemOnManifold::BoxTrajProblemOnManifold(
    const Manifold& M, const std::string& configPath)
    : Problem(M), config_(configPath)
{
  nBoxes_ = static_cast<size_t>(config_["nBoxes"].asInt());
  nObstacles_ = static_cast<size_t>(config_["nObstacles"].asInt());
  nPlans_ = nBoxes_ * nObstacles_;

  boxSize_ = config_["BoxSize"].asVector3d();
  initPos_ = config_["initPos"].asVector3d();
  finalPos_ = config_["finalPos"].asVector3d();

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

  for (size_t i = 0; i < nBoxes_; ++i)
  {
    for (size_t j = 0; j < nObstacles_; ++j)
    {
      plans_.push_back(Plan(j, i));
    }
  }

  for (size_t i = 0; i < nPlans_; i++)
  {
    cstrNames_.push_back("Plan" + std::to_string(i) + "BetweenBox" +
                         std::to_string(plans_[i].boxAbove()) + "AndObstacle" +
                         std::to_string(plans_[i].boxBelow()));
  }

  std::stringstream sstm;
  sstm << "FeetTraj" << nBoxes_ << "Boxes";
  name() = sstm.str();

  outRepObjDiff_.resize(1, M.representationDim());
  outRepObjDiff_.setZero();
  outRep_.resize(static_cast<Index>(16 * nPlans_), M.representationDim());
  // outRep_.resize(
  // static_cast<Index>(8 * boxAboveFixedPlanCstrs_.size() + 16 * nPlans_),
  // M.representationDim());
  outRep_.setZero();
}

BoxTrajProblemOnManifold::~BoxTrajProblemOnManifold() {}

Eigen::VectorXd BoxTrajProblemOnManifold::findInitPoint()
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

CartesianProduct* BoxTrajProblemOnManifold::buildManifold(
    const Index& nBoxes, const Index& nObstacles)
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

  Index nPlanes = nBoxes * nObstacles;
  CartesianProduct* MPlanes;
  if (nPlanes == 0)
    MPlanes = new CartesianProduct(*R0, *R0);
  else if (nPlanes == 1)
    MPlanes = new CartesianProduct(*MPlane, *R0);
  else
  {
    MPlanes = new CartesianProduct(*MPlane, *MPlane);
    for (size_t i = 2; i < static_cast<size_t>(nPlanes); ++i) MPlanes->multiply(*MPlane);
  }
  CartesianProduct* MBoxesAndPlanes = new CartesianProduct(*MBoxes, *MPlanes);
  return MBoxesAndPlanes;
}

void BoxTrajProblemOnManifold::getTangentLB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = -infinity;
}

void BoxTrajProblemOnManifold::getTangentUB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = infinity;
}

void BoxTrajProblemOnManifold::evalObj(double& out) const
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
void BoxTrajProblemOnManifold::evalObjDiff(RefMat out) const
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
  outRepObjDiff_.block(0, (static_cast<Index>(nBoxes_) - 1) * boxRepDim, 1, boxRepDim) +=
      -2 * dist.transpose();
  M().applyDiffRetractation(out, outRepObjDiff_, x().value());
}

void BoxTrajProblemOnManifold::evalLinCstr(RefVec, size_t) const
{
  // if (i == 0)  // Initial point constraint
  // out << phi_x_z()(0)(0)[0];
  // else if (i == 1)  // Final point constraint
  // out << phi_x_z()(0)(static_cast<Index>(nBoxes_) - 1)[0];
}

void BoxTrajProblemOnManifold::evalLinCstrDiff(RefMat, size_t) const
{
  // if (i == 0)  // Initial point constraint
  // out.block(0, 0, 3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  // else if (i == 1)  // Final point constraint
  // out.block(0, 3 * (nBoxes_ - 1), 3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1;
}

void BoxTrajProblemOnManifold::getLinCstrLB(RefVec, size_t) const
{
  // if (i == 0)  // Initial point constraint
  // out << initPos_;
  // else if (i == 1)  // Final point constraint
  // out << finalPos_;
}

void BoxTrajProblemOnManifold::getLinCstrUB(RefVec, size_t) const
{
  // if (i == 0)  // Initial point constraint
  // out << initPos_;
  // else if (i == 1)  // Final point constraint
  // out << finalPos_;
}

void BoxTrajProblemOnManifold::evalNonLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  out.setZero();

  // if (i < 2)
  //{
  //// nothing to do
  //}
  // else
  //{
  auto iPlan = i;
  auto iBoxAbove = plans_[iPlan].boxAbove();
  auto iBoxBelow = plans_[iPlan].boxBelow();
  Eigen::Vector3d transAbove = phi_x_z()(0)(iBoxAbove)[0];
  Eigen::Vector4d quatAbove(0, 0, 0, 1);
  Eigen::Vector3d transBelow = obstacles_[iBoxBelow].center();
  Eigen::Vector4d quatBelow(0, 0, 0, 1);
  Eigen::Vector3d normal = phi_x_z()(1)(iPlan)[1];
  double d = phi_x_z()(1)(iPlan)[0](0);
  std::cout << "\nCompute cstr box above: " << std::endl;
  boxAbovePlanFcts_[iBoxAbove].compute(out.head(8), transAbove, quatAbove, d,
                                       normal);
  std::cout << "\nCompute cstr obstacle below: " << std::endl;
  obstacleAbovePlanFcts_[iBoxBelow].compute(out.tail(8), transBelow, quatBelow,
                                            -d, -normal);
  //}
}

void BoxTrajProblemOnManifold::evalNonLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == nonLinCstrDim(i) && "wrong row size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");

  // if (i < 2)
  //{
  //// nothing to do
  //}
  // else
  //{
  auto iPlan = i;
  auto iBoxAbove = plans_[iPlan].boxAbove();
  auto iBoxBelow = plans_[iPlan].boxBelow();
  Eigen::Vector3d transAbove = phi_x_z()(0)(iBoxAbove)[0];
  Eigen::Vector4d quatAbove(1, 0, 0, 0);
  Eigen::Vector3d transBelow = obstacles_[iBoxBelow].center();
  Eigen::Vector4d quatBelow(1, 0, 0, 0);
  Eigen::Vector3d normal = phi_x_z()(1)(iPlan)[1];
  double d = phi_x_z()(1)(iPlan)[0](0);
  auto rowBegin = 16 * iPlan;
  auto boxAboveBegin = 3 * iBoxAbove;
  auto planBegin = 3 * nBoxes_ + 4 * iPlan;
  long rowBeginLong = static_cast<long>(rowBegin);
  long boxAboveBeginLong = static_cast<long>(boxAboveBegin);
  long planBeginLong = static_cast<long>(planBegin);

  Eigen::Matrix<double, 8, 1> tmpWTF;

  boxAbovePlanFcts_[iBoxAbove].diffTrans(
      outRep_.block(rowBeginLong, boxAboveBeginLong, 8, 3), transAbove,
      quatAbove, d, normal);
  boxAbovePlanFcts_[iBoxAbove].diffD(tmpWTF, transAbove, quatAbove, d, normal);
  outRep_.block(rowBeginLong, planBeginLong, 8, 1) = tmpWTF;
  boxAbovePlanFcts_[iBoxAbove].diffNormal(
      outRep_.block(rowBeginLong, planBeginLong + 1, 8, 3), transAbove,
      quatAbove, d, normal);

  obstacleAbovePlanFcts_[iBoxBelow].diffD(tmpWTF, transBelow, quatBelow, -d,
                                          -normal);
  outRep_.block(rowBeginLong + 8, planBeginLong, 8, 1) = -tmpWTF;

  Eigen::Matrix<double, 8, 3> tmpDiffNormal;
  obstacleAbovePlanFcts_[iBoxBelow].diffNormal(tmpDiffNormal, transBelow,
                                               quatBelow, -d, -normal);
  outRep_.block(rowBeginLong + 8, planBeginLong + 1, 8, 3) = -tmpDiffNormal;

  M().applyDiffRetractation(out, outRep_.middleRows(rowBeginLong, 16),
                            x().value());
  //}
}

void BoxTrajProblemOnManifold::getNonLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  // if (i < 2)
  //{
  //// nothing to do
  //}
  // else
  //{
  auto iPlan = i;
  auto iBoxAbove = plans_[iPlan].boxAbove();
  auto iBoxBelow = plans_[iPlan].boxBelow();
  boxAbovePlanFcts_[iBoxAbove].LB(out.head(8));
  boxAbovePlanFcts_[iBoxBelow].LB(out.tail(8));
  //}
}
void BoxTrajProblemOnManifold::getNonLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  // if (i < 2)
  //{
  //// nothing to do
  //}
  // else
  //{
  auto iPlan = i;
  auto iBoxAbove = plans_[iPlan].boxAbove();
  auto iBoxBelow = plans_[iPlan].boxBelow();
  boxAbovePlanFcts_[iBoxAbove].UB(out.head(8));
  boxAbovePlanFcts_[iBoxBelow].UB(out.tail(8));
  //}
}

size_t BoxTrajProblemOnManifold::numberOfCstr() const { return nPlans_; }

Index BoxTrajProblemOnManifold::linCstrDim(size_t) const
{
  // if (i < 2)
  // return 3;
  // else
  return 0;
}

Index BoxTrajProblemOnManifold::nonLinCstrDim(size_t) const
{
  // if (i < 2)
  // return 0;
  // else
  return 16;
}

std::string BoxTrajProblemOnManifold::getCstrName(const size_t) const
{
  std::string str("");
  return str;
}

// void BoxTrajProblemOnManifold::fileForMatlab(std::string fileName,
// const mnf::Point& x) const
//{
//}
} /* feettrajectory */
