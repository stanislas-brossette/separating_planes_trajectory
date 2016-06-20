#include <iostream>
#include <fstream>
#include <limits>

#include <Eigen/Geometry>

#include <cube-stacks/CubeStackProblemOnManifold.hh>
#include <cube-stacks/utils/ProblemConfig.hh>
#include <cube-stacks/utils/quat2mat.hh>

#include <pgsolver/utils/usingManifold.h>

#include <manifolds/RealSpace.h>
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>
using mnf::RealSpace;
using mnf::SO3;
using mnf::S2;
using mnf::CartesianProduct;
using mnf::ExpMapMatrix;
using mnf::ExpMapQuaternion;
namespace cubestacks
{
CubeStackProblemOnManifold::CubeStackProblemOnManifold(
    const Manifold& M, const std::string& configPath)
    : Problem(M), config_(configPath)
{
  nCubes_ = static_cast<size_t>(M(0).dim()) / 6;

  normalZPlus_ = config_["normalZPlus"];
  normalXMinus_ = config_["normalXMinus"];
  normalXPlus_ = config_["normalXPlus"];
  normalYMinus_ = config_["normalYMinus"];
  normalYPlus_ = config_["normalYPlus"];

  distZPlus_ = config_["distZPlus"];
  distXMinus_ = config_["distXMinus"];
  distXPlus_ = config_["distXPlus"];
  distYMinus_ = config_["distYMinus"];
  distYPlus_ = config_["distYPlus"];

  Eigen::VectorXd lCubes = config_["lCubes"].asVectorXd();

  for (int i = 0; i < static_cast<int>(nCubes_); ++i)
  {
    Cube c(i, lCubes(i));
    cubeAbovePlanFcts_.push_back(CubeAbovePlan(c));
    cubes_.push_back(c);
    cubeAboveFixedPlanCstrs_.push_back(
        CubeAboveFixedPlan(c, normalZPlus_, distZPlus_));
    cstrNames_.push_back("Cube" + std::to_string(i) + "AboveFixedPlan" +
                         "ZPlus");
    cubeAboveFixedPlanCstrs_.push_back(
        CubeAboveFixedPlan(c, normalXPlus_, distXPlus_));
    cstrNames_.push_back("Cube" + std::to_string(i) + "AboveFixedPlan" +
                         "XPlus");
    cubeAboveFixedPlanCstrs_.push_back(
        CubeAboveFixedPlan(c, normalXMinus_, distXMinus_));
    cstrNames_.push_back("Cube" + std::to_string(i) + "AboveFixedPlan" +
                         "XMinus");
    cubeAboveFixedPlanCstrs_.push_back(
        CubeAboveFixedPlan(c, normalYPlus_, distYPlus_));
    cstrNames_.push_back("Cube" + std::to_string(i) + "AboveFixedPlan" +
                         "YPlus");
    cubeAboveFixedPlanCstrs_.push_back(
        CubeAboveFixedPlan(c, normalYMinus_, distYMinus_));
    cstrNames_.push_back("Cube" + std::to_string(i) + "AboveFixedPlan" +
                         "YMinus");
  }

  nPlans_ = (nCubes_ * (nCubes_ - 1)) / 2;

  for (size_t i = 0; i < nCubes_; ++i)
  {
    for (size_t j = i + 1; j < nCubes_; ++j)
    {
      plans_.push_back(Plan(i, j));
    }
  }

  for (size_t i = 0; i < nPlans_; i++)
  {
    cstrNames_.push_back("Plan" + std::to_string(i) + "BetweenCube" +
                         std::to_string(plans_[i].cubeAbove()) + "AndCube" +
                         std::to_string(plans_[i].cubeBelow()));
  }

  std::stringstream sstm;
  sstm << "Stack" << nCubes_ << "CubesS2";
  name() = sstm.str();

  outRepObjDiff_.resize(1, M.representationDim());
  outRepObjDiff_.setZero();
  outRep_.resize(static_cast<Index>(8 * cubeAboveFixedPlanCstrs_.size() + 16*nPlans_),
                 M.representationDim());
  outRep_.setZero();
}

Eigen::VectorXd CubeStackProblemOnManifold::findInitPoint()
{
  Point xRand = M().createRandomPoint();
  Eigen::Vector3d interCubes;
  for (size_t i = 0; i < nCubes_; ++i)
  {
    xRand(0)(i)[0][2] += 2 * static_cast<double>(i);
  }
  for (size_t i = 0; i < nPlans_; ++i)
  {
    size_t cubeBelow = plans_[i].cubeBelow();
    size_t cubeAbove = plans_[i].cubeAbove();
    Eigen::Vector3d cubeBelowPos = xRand(0)(cubeBelow)[0];
    Eigen::Vector3d cubeAbovePos = xRand(0)(cubeAbove)[0];
    interCubes = cubeAbovePos - cubeBelowPos;
    interCubes = interCubes / interCubes.norm();
    xRand(1)(i)[1] = interCubes;
    xRand(1)(i)[0][0] = interCubes.dot((cubeAbovePos + cubeBelowPos) / 2);
  }
  return xRand.value();
}

CartesianProduct* CubeStackProblemOnManifold::buildManifold(Index nCubes)
{
  assert(nCubes > 0 && "Number of Cubes must be positive and non null");
  RealSpace* R1_ = new RealSpace(1);
  RealSpace* R3_ = new RealSpace(3);
  SO3<ExpMapQuaternion>* SO3_ = new SO3<ExpMapQuaternion>();
  S2* S2_ = new S2();
  CartesianProduct* MCube = new CartesianProduct(*R3_, *SO3_);
  CartesianProduct* MPlane = new CartesianProduct(*R1_, *S2_);
  RealSpace* R0 = new RealSpace(0);
  CartesianProduct* MCubes;
  if (nCubes == 1)
    MCubes = new CartesianProduct(*MCube, *R0);
  else
  {
    MCubes = new CartesianProduct(*MCube, *MCube);
    for (size_t i = 2; i < static_cast<size_t>(nCubes); ++i)
      MCubes->multiply(*MCube);
  }

  Index nPlanes = nCubes * (nCubes - 1) / 2;
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
  CartesianProduct* MCubesAndPlanes = new CartesianProduct(*MCubes, *MPlanes);
  return MCubesAndPlanes;
}

void CubeStackProblemOnManifold::getTangentLB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = -infinity;
}
void CubeStackProblemOnManifold::getTangentUB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = infinity;
}

void CubeStackProblemOnManifold::evalObj(double& out) const
{
  out = 0;
  for (size_t i = 0; i < nCubes_; ++i)
  {
    Eigen::Vector3d pos = phi_x_z()(0)(i)[0];
    out += pos[2];
  }
}
void CubeStackProblemOnManifold::evalObjDiff(RefMat out) const
{
  assert(out.rows() == 1 && "wrong rows size");
  assert(out.cols() == M().tangentDim() && "wrong cols size");
  outRepObjDiff_.setZero();
  for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
  {
    Index cubeRepDim = M()(0)(static_cast<size_t>(i)).representationDim();
    outRepObjDiff_.block(0, i * cubeRepDim, 1, cubeRepDim) << 0, 0, 1, 0, 0, 0,
        0;
  }
  M().applyDiffRetractation(out, outRepObjDiff_, x().value());
}

void CubeStackProblemOnManifold::evalLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
}

void CubeStackProblemOnManifold::evalLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == linCstrDim(i) && "Wrong rows size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");
}

void CubeStackProblemOnManifold::getLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
}

void CubeStackProblemOnManifold::getLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
}

void CubeStackProblemOnManifold::evalNonLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  out.setZero();

  if (i < cubeAboveFixedPlanCstrs_.size())
  {
    unsigned long iC =
        static_cast<unsigned long>(cubeAboveFixedPlanCstrs_[i].cube().index());
    Eigen::Vector3d trans = phi_x_z()(0)(iC)[0];
    Eigen::Vector4d quat = phi_x_z()(0)(iC)[1];
    cubeAboveFixedPlanCstrs_[i].compute(out, trans, quat);
  }
  else
  {
    auto iPlan = i - cubeAboveFixedPlanCstrs_.size();
    auto iCubeAbove = plans_[iPlan].cubeAbove();
    auto iCubeBelow = plans_[iPlan].cubeBelow();
    Eigen::Vector3d transAbove = phi_x_z()(0)(iCubeAbove)[0];
    Eigen::Vector4d quatAbove = phi_x_z()(0)(iCubeAbove)[1];
    Eigen::Vector3d transBelow = phi_x_z()(0)(iCubeBelow)[0];
    Eigen::Vector4d quatBelow = phi_x_z()(0)(iCubeBelow)[1];
    Eigen::Vector3d normal = phi_x_z()(1)(iPlan)[1];
    double d = phi_x_z()(1)(iPlan)[0](0);
    cubeAbovePlanFcts_[iCubeAbove].compute(out.head(8), transAbove, quatAbove,
                                           d, normal);
    cubeAbovePlanFcts_[iCubeBelow].compute(out.tail(8), transBelow, quatBelow,
                                           -d, -normal);
  }
}

void CubeStackProblemOnManifold::evalNonLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == nonLinCstrDim(i) && "wrong row size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");

  if (i < cubeAboveFixedPlanCstrs_.size())
  {
    int iC = cubeAboveFixedPlanCstrs_[i].cube().index();
    Eigen::Vector3d trans = phi_x_z()(0)(static_cast<size_t>(iC))[0];
    Eigen::Vector4d quat = phi_x_z()(0)(static_cast<size_t>(iC))[1];
    cubeAboveFixedPlanCstrs_[i].diffTrans(
        outRep_.block(8 * static_cast<long>(i), 7 * iC, 8, 3), trans, quat);
    cubeAboveFixedPlanCstrs_[i].diffQuat(
        outRep_.block(8 * static_cast<long>(i), 7 * iC + 3, 8, 4), trans, quat);
    M().applyDiffRetractation(
        out, outRep_.middleRows(8 * static_cast<long>(i), 8), x().value());
  }
  else
  {
    auto iPlan = i - cubeAboveFixedPlanCstrs_.size();
    auto iCubeAbove = plans_[iPlan].cubeAbove();
    auto iCubeBelow = plans_[iPlan].cubeBelow();
    Eigen::Vector3d transAbove = phi_x_z()(0)(iCubeAbove)[0];
    Eigen::Vector4d quatAbove = phi_x_z()(0)(iCubeAbove)[1];
    Eigen::Vector3d transBelow = phi_x_z()(0)(iCubeBelow)[0];
    Eigen::Vector4d quatBelow = phi_x_z()(0)(iCubeBelow)[1];
    Eigen::Vector3d normal = phi_x_z()(1)(iPlan)[1];
    double d = phi_x_z()(1)(iPlan)[0](0);
    auto rowBegin = 8 * cubeAboveFixedPlanCstrs_.size() + 16 * iPlan;
    auto cubeAboveBegin = 7 * iCubeAbove;
    auto cubeBelowBegin = 7 * iCubeBelow;
    auto planBegin = 7 * nCubes_ + 4 * iPlan;
    long rowBeginLong = static_cast<long>(rowBegin);
    long cubeAboveBeginLong = static_cast<long>(cubeAboveBegin);
    long cubeBelowBeginLong = static_cast<long>(cubeBelowBegin);
    long planBeginLong = static_cast<long>(planBegin);

    Eigen::Matrix<double, 8, 1> tmpWTF;

    cubeAbovePlanFcts_[iCubeAbove].diffTrans(
        outRep_.block(rowBeginLong, cubeAboveBeginLong, 8, 3), transAbove,
        quatAbove, d, normal);
    cubeAbovePlanFcts_[iCubeAbove].diffQuat(
        outRep_.block(rowBeginLong, cubeAboveBeginLong + 3, 8, 4), transAbove,
        quatAbove, d, normal);
    cubeAbovePlanFcts_[iCubeAbove].diffD(tmpWTF, transAbove, quatAbove, d,
                                         normal);
    outRep_.block(rowBeginLong, planBeginLong, 8, 1) = tmpWTF;
    cubeAbovePlanFcts_[iCubeAbove].diffNormal(
        outRep_.block(rowBeginLong, planBeginLong + 1, 8, 3), transAbove,
        quatAbove, d, normal);

    cubeAbovePlanFcts_[iCubeBelow].diffTrans(
        outRep_.block(rowBeginLong + 8, cubeBelowBeginLong, 8, 3), transBelow,
        quatBelow, -d, -normal);
    cubeAbovePlanFcts_[iCubeBelow].diffQuat(
        outRep_.block(rowBeginLong + 8, cubeBelowBeginLong + 3, 8, 4),
        transBelow, quatBelow, -d, -normal);
    cubeAbovePlanFcts_[iCubeBelow].diffD(tmpWTF, transBelow, quatBelow, -d,
                                         -normal);
    outRep_.block(rowBeginLong + 8, planBeginLong, 8, 1) = -tmpWTF;

    Eigen::Matrix<double, 8, 3> tmpDiffNormal;
    cubeAbovePlanFcts_[iCubeBelow].diffNormal(
        tmpDiffNormal, transBelow,
        quatBelow, -d, -normal);
    outRep_.block(rowBeginLong + 8, planBeginLong + 1, 8, 3) = -tmpDiffNormal;

    M().applyDiffRetractation(out, outRep_.middleRows(rowBeginLong, 16),
                              x().value());
  }
}

void CubeStackProblemOnManifold::getNonLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  if (i < cubeAboveFixedPlanCstrs_.size())
  {
    cubeAboveFixedPlanCstrs_[i].LB(out);
  }
  else
  {
    auto iPlan = i - cubeAboveFixedPlanCstrs_.size();
    auto iCubeAbove = plans_[iPlan].cubeAbove();
    auto iCubeBelow = plans_[iPlan].cubeBelow();
    cubeAbovePlanFcts_[iCubeAbove].LB(out.head(8));
    cubeAbovePlanFcts_[iCubeBelow].LB(out.tail(8));
  }
}
void CubeStackProblemOnManifold::getNonLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  if (i < cubeAboveFixedPlanCstrs_.size())
  {
    cubeAboveFixedPlanCstrs_[i].UB(out);
  }
  else
  {
    auto iPlan = i - cubeAboveFixedPlanCstrs_.size();
    auto iCubeAbove = plans_[iPlan].cubeAbove();
    auto iCubeBelow = plans_[iPlan].cubeBelow();
    cubeAbovePlanFcts_[iCubeAbove].UB(out.head(8));
    cubeAbovePlanFcts_[iCubeBelow].UB(out.tail(8));
  }
}

size_t CubeStackProblemOnManifold::numberOfCstr() const
{
  size_t nb = cubeAboveFixedPlanCstrs_.size();
  nb += nPlans_;
  return nb;
}

Index CubeStackProblemOnManifold::linCstrDim(size_t) const { return 0; }

Index CubeStackProblemOnManifold::nonLinCstrDim(size_t i) const
{
  if (i < cubeAboveFixedPlanCstrs_.size())
    return 8;
  else
    return 16;
}

std::string CubeStackProblemOnManifold::getCstrName(const size_t i) const
{
  return cstrNames_[i];
}

void CubeStackProblemOnManifold::fileForMatlab(std::string fileName,
                                       const Point& x) const
{
  std::ofstream myFile;
  myFile.open(fileName);

  myFile << "clear all;" << std::endl;
  myFile << "clf;" << std::endl;
  myFile << "material shiny;" << std::endl;
  myFile << "alpha(\'color\');" << std::endl;
  myFile << "alphamap(\'rampdown\');" << std::endl;
  myFile << "xlabel(\'x\');" << std::endl;
  myFile << "ylabel(\'y\');" << std::endl;
  myFile << "zlabel(\'z\');" << std::endl;
  myFile << "grid \'on\';" << std::endl;
  myFile << "view(30,30);" << std::endl;
  myFile << "axis(\'equal\');   " << std::endl;

  double minX(-2.0);
  double maxX(2.0);
  double minY(-2.0);
  double maxY(2.0);
  double minZ(-2.0);
  double maxZ(2.0);
  for (size_t i = 0; i < nCubes_; ++i)
  {
    Eigen::Vector3d pC = x(0)(i)[0];
    minX = fmin(minX, pC(0) - cubes_[i].l());
    maxX = fmax(maxX, pC(0) + cubes_[i].l());
    minY = fmin(minY, pC(1) - cubes_[i].l());
    maxY = fmax(maxY, pC(1) + cubes_[i].l());
    minZ = fmin(minZ, pC(2) - cubes_[i].l());
    maxZ = fmax(maxZ, pC(2) + cubes_[i].l());
    Eigen::Matrix3d rotC = quat2mat(x(0)(i)[1]);
    myFile << "vert" << i << " = [ ";
    for (Index j = 0; j < 8; ++j)
    {
      Eigen::Vector3d v =
          pC + rotC * (cubes_[i].vertex(static_cast<size_t>(j)));
      myFile << v.transpose();

      if (j < 7)
      {
        myFile << "; ..." << std::endl;
      }
      else
      {
        myFile << "];" << std::endl;
      }
    }
    myFile << "fac" << i
           << " = [1 2 4 3; 2 6 8 4; 4 3 7 8; 1 5 7 3; 1 2 6 5; 5 6 8 7];"
           << std::endl;
    myFile << "patch( \'Faces\', fac" << i << ", \'Vertices\', vert" << i
           << ", \'FaceColor\', \'r\');" << std::endl;
  }
  for (size_t i = 0; i < nPlans_; ++i)
  {
    double r = x(1)(i)[0][0];
    Eigen::Vector3d n = x(1)(i)[1];
    size_t cubeBelow = plans_[i].cubeBelow();
    size_t cubeAbove = plans_[i].cubeAbove();
    Eigen::Vector3d Cb = x(0)(cubeBelow)[0];
    Eigen::Vector3d Cu = x(0)(cubeAbove)[0];

    double lambda = (r - n.dot(Cb)) / n.dot(Cu - Cb);
    Eigen::Vector3d zeroPlane = Cb + lambda * (Cu - Cb);

    Eigen::Vector3d t(n(1), -n(0), 0);
    t = t / t.norm();

    Eigen::Vector3d b = n.cross(t);
    b = b / b.norm();

    myFile << "vert" << i << " = [ ";

    Eigen::Vector3d vertex;
    vertex = zeroPlane + t + b;
    myFile << vertex.transpose() << "; ..." << std::endl;
    vertex = zeroPlane + t - b;
    myFile << vertex.transpose() << "; ..." << std::endl;
    vertex = zeroPlane - t - b;
    myFile << vertex.transpose() << "; ..." << std::endl;
    vertex = zeroPlane - t + b;
    myFile << vertex.transpose() << "; ..." << std::endl;
    myFile << "];" << std::endl;
    myFile << "fac" << i << " = [1 2 3 4];" << std::endl;
    myFile << "patch( \'Faces\', fac" << i << ", \'Vertices\', vert" << i
           << ", \'FaceColor\', \'y\', 'FaceAlpha', 0.5);" << std::endl;
  }
  myFile << "xlim([" << minX << " " << maxX << "]);" << std::endl;
  myFile << "ylim([" << minY << " " << maxY << "]);" << std::endl;
  myFile << "zlim([" << minZ << " " << maxZ << "]);" << std::endl;

  myFile << "vertPlanZPlus = [ ";
  myFile << -distXMinus_ << ", " << -distYMinus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << -distXMinus_ << ", " << distYPlus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << distYPlus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << -distYMinus_ << ", " << distZPlus_ << "];"
         << std::endl;
  myFile << "facZPlus = [1 2 3 4];" << std::endl;
  myFile << "patch( \'Faces\', facZPlus, \'Vertices\', vertPlanZPlus, "
            "\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  myFile << "vertPlanXPlus = [ ";
  myFile << distXPlus_ << ", " << -distYMinus_ << ", " << distYPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << distYPlus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << distYPlus_ << ", " << distZPlus_ + 3 << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << -distYMinus_ << ", " << distZPlus_ + 3 << "];"
         << std::endl;
  myFile << "facXPlus = [1 2 3 4];" << std::endl;
  myFile << "patch( \'Faces\', facXPlus, \'Vertices\', vertPlanXPlus, "
            "\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  myFile << "vertPlanXMinus = [ ";
  myFile << -distXMinus_ << ", " << -distYMinus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << -distXMinus_ << ", " << distYPlus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << -distXMinus_ << ", " << distYPlus_ << ", " << distZPlus_ + 3
         << ";..." << std::endl;
  myFile << -distXMinus_ << ", " << -distYMinus_ << ", " << distZPlus_ + 3
         << "];" << std::endl;
  myFile << "facXMinus = [1 2 3 4];" << std::endl;
  myFile << "patch( \'Faces\', facXMinus, \'Vertices\', vertPlanXMinus, "
            "\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  myFile << "vertPlanYPlus = [ ";
  myFile << -distXMinus_ << ", " << distYPlus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << distYPlus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << distYPlus_ << ", " << distZPlus_ + 3 << ";..."
         << std::endl;
  myFile << -distXMinus_ << ", " << distYPlus_ << ", " << distZPlus_ + 3 << "];"
         << std::endl;
  myFile << "facYPlus = [1 2 3 4];" << std::endl;
  myFile << "patch( \'Faces\', facYPlus, \'Vertices\', vertPlanYPlus, "
            "\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  myFile << "vertPlanYMinus = [ ";
  myFile << -distXMinus_ << ", " << -distYMinus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << -distYMinus_ << ", " << distZPlus_ << ";..."
         << std::endl;
  myFile << distXPlus_ << ", " << -distYMinus_ << ", " << distZPlus_ + 3
         << ";..." << std::endl;
  myFile << -distXMinus_ << ", " << -distYMinus_ << ", " << distZPlus_ + 3
         << "];" << std::endl;
  myFile << "facYMinus = [1 2 3 4];" << std::endl;
  myFile << "patch( \'Faces\', facYMinus, \'Vertices\', vertPlanYMinus, "
            "\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  myFile.close();
}
}
