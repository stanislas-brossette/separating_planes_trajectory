#include <iostream>
#include <fstream>
#include <limits>

#include <Eigen/Geometry>

#include <cube-stacks/CubeStackProblemOnR.hh>
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
CubeStackProblemOnR::CubeStackProblemOnR(
    const Manifold& M, const std::string& configPath)
    : Problem(M), config_(configPath)
{
  nCubes_ = static_cast<size_t>(config_["nCubes"].asInt());

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

  std::cout << "nCubes_ = \n" << nCubes_ << std::endl;
  std::cout << "nPlans_ = \n" << nPlans_ << std::endl;

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
  for (size_t i = 0; i < nCubes_; i++) 
  {
    cstrNames_.push_back("Norm1QuatCube" + std::to_string(i));
  }
  for (size_t i = 0; i < nPlans_; i++) 
  {
    cstrNames_.push_back("Norm1S2Plan" + std::to_string(i));
  }

  std::stringstream sstm;
  sstm << "Stack" << nCubes_ << "Cubes_NO_Manifold";
  name() = sstm.str();
  
  iMaxCstrFixPlan_ = cubeAboveFixedPlanCstrs_.size();
  iMaxCstrS2Plan_ = iMaxCstrFixPlan_ + nPlans_;
  iMaxCstrQuatNorm1_ = iMaxCstrS2Plan_ + nCubes_;
  iMaxCstrVecNorm1_ = iMaxCstrQuatNorm1_ + nPlans_;
}

Eigen::VectorXd CubeStackProblemOnR::findInitPoint()
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

CartesianProduct* CubeStackProblemOnR::buildManifold(Index nCubes)
{
  assert(nCubes > 0 && "Number of Cubes must be positive and non null");
  RealSpace* R1_ = new RealSpace(1);
  RealSpace* R3_ = new RealSpace(3);
  RealSpace* SO3_ = new RealSpace(4);
  RealSpace* S2_ = new RealSpace(3);
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

void CubeStackProblemOnR::getTangentLB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = -infinity;
}
void CubeStackProblemOnR::getTangentUB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) out(i) = infinity;
}

void CubeStackProblemOnR::evalObj(double& out) const
{
  out = 0;
  for (size_t i = 0; i < nCubes_; ++i)
  {
    Eigen::Vector3d pos = phi_x_z()(0)(i)[0];
    out += pos[2];
  }
}
void CubeStackProblemOnR::evalObjDiff(RefMat out) const
{
  assert(out.rows() == 1 && "wrong rows size");
  assert(out.cols() == M().tangentDim() && "wrong cols size");
  for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
  {
    Index cubeRepDim = M()(0)(static_cast<size_t>(i)).representationDim();
    out.block(0, i * cubeRepDim, 1, cubeRepDim) << 0, 0, 1, 0, 0, 0,
        0;
  }
}

void CubeStackProblemOnR::evalLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
}

void CubeStackProblemOnR::evalLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == linCstrDim(i) && "Wrong rows size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");
}

void CubeStackProblemOnR::getLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
}

void CubeStackProblemOnR::getLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
}

void CubeStackProblemOnR::evalNonLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  out.setZero();

  if (i < iMaxCstrFixPlan_)
  {
    unsigned long iC =
        static_cast<unsigned long>(cubeAboveFixedPlanCstrs_[i].cube().index());
    Eigen::Vector3d trans = phi_x_z()(0)(iC)[0];
    Eigen::Vector4d quat = phi_x_z()(0)(iC)[1];
    cubeAboveFixedPlanCstrs_[i].compute(out, trans, quat);
  }
  else if(iMaxCstrFixPlan_ <= i && i < iMaxCstrS2Plan_)
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
  else if(iMaxCstrS2Plan_ <= i && i < iMaxCstrQuatNorm1_)
  {
    auto iCube = i - iMaxCstrS2Plan_;
    Eigen::Vector4d quat = phi_x_z()(0)(iCube)[1];
    out(0) = quat.squaredNorm();
  }
  else if(iMaxCstrQuatNorm1_ <= i && i < iMaxCstrVecNorm1_)
  {
    auto iPlan = i - iMaxCstrQuatNorm1_;
    Eigen::Vector3d vec = phi_x_z()(1)(iPlan)[1];
    out(0) = vec.squaredNorm();
  }
  else
  {
    std::cerr << "WOOOPS! Wrong csrt index " << i << std::endl;
  }
}

void CubeStackProblemOnR::evalNonLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == nonLinCstrDim(i) && "wrong row size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");

  if (i < iMaxCstrFixPlan_)
  {
    int iC = cubeAboveFixedPlanCstrs_[i].cube().index();
    Eigen::Vector3d trans = phi_x_z()(0)(static_cast<size_t>(iC))[0];
    Eigen::Vector4d quat = phi_x_z()(0)(static_cast<size_t>(iC))[1];
    cubeAboveFixedPlanCstrs_[i].diffTrans(
        out.block(0, 7 * iC, 8, 3), trans, quat);
    cubeAboveFixedPlanCstrs_[i].diffQuat(
        out.block(0, 7 * iC + 3, 8, 4), trans, quat);
  }
  else if(iMaxCstrFixPlan_ <= i && i < iMaxCstrS2Plan_)
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
    auto rowBegin = 0; //8 * cubeAboveFixedPlanCstrs_.size() + 16 * iPlan;
    auto cubeAboveBegin = 7 * iCubeAbove;
    auto cubeBelowBegin = 7 * iCubeBelow;
    auto planBegin = 7 * nCubes_ + 4 * iPlan;
    long rowBeginLong = static_cast<long>(rowBegin);
    long cubeAboveBeginLong = static_cast<long>(cubeAboveBegin);
    long cubeBelowBeginLong = static_cast<long>(cubeBelowBegin);
    long planBeginLong = static_cast<long>(planBegin);

    Eigen::Matrix<double, 8, 1> tmpWTF;

    cubeAbovePlanFcts_[iCubeAbove].diffTrans(
        out.block(rowBeginLong, cubeAboveBeginLong, 8, 3), transAbove,
        quatAbove, d, normal);
    cubeAbovePlanFcts_[iCubeAbove].diffQuat(
        out.block(rowBeginLong, cubeAboveBeginLong + 3, 8, 4), transAbove,
        quatAbove, d, normal);
    cubeAbovePlanFcts_[iCubeAbove].diffD(tmpWTF, transAbove, quatAbove, d,
                                         normal);
    out.block(rowBeginLong, planBeginLong, 8, 1) = tmpWTF;
    cubeAbovePlanFcts_[iCubeAbove].diffNormal(
        out.block(rowBeginLong, planBeginLong + 1, 8, 3), transAbove,
        quatAbove, d, normal);

    cubeAbovePlanFcts_[iCubeBelow].diffTrans(
        out.block(rowBeginLong + 8, cubeBelowBeginLong, 8, 3), transBelow,
        quatBelow, -d, -normal);
    cubeAbovePlanFcts_[iCubeBelow].diffQuat(
        out.block(rowBeginLong + 8, cubeBelowBeginLong + 3, 8, 4),
        transBelow, quatBelow, -d, -normal);
    cubeAbovePlanFcts_[iCubeBelow].diffD(tmpWTF, transBelow, quatBelow, -d,
                                         -normal);
    out.block(rowBeginLong + 8, planBeginLong, 8, 1) = -tmpWTF;

    Eigen::Matrix<double, 8, 3> tmpDiffNormal;
    cubeAbovePlanFcts_[iCubeBelow].diffNormal(
        tmpDiffNormal, transBelow,
        quatBelow, -d, -normal);
    out.block(rowBeginLong + 8, planBeginLong + 1, 8, 3) = -tmpDiffNormal;
  }
  else if(iMaxCstrS2Plan_ <= i && i < iMaxCstrQuatNorm1_)
  {
    auto iCube = i - iMaxCstrS2Plan_;
    auto quatBegin = 7 * iCube + 3;
    Eigen::Vector4d q = phi_x_z()(0)(iCube)[1];
    out.middleCols(static_cast<long>(quatBegin), 4) << 2*q(0), 2*q(1), 2*q(2), 2*q(3);
  }
  else if(iMaxCstrQuatNorm1_ <= i && i < iMaxCstrVecNorm1_)
  {
    auto iPlan = i - iMaxCstrQuatNorm1_;
    auto vecBegin = 7 * nCubes_ + 4 * iPlan + 1;
    Eigen::Vector3d v = phi_x_z()(1)(iPlan)[1];
    out.middleCols(static_cast<long>(vecBegin), 3) << 2*v(0), 2*v(1), 2*v(2);
  }
  else
  {
    std::cerr << "WOOOPS! Wrong csrt index " << i << std::endl;
  }
}

void CubeStackProblemOnR::getNonLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  if (i < iMaxCstrFixPlan_)
  {
    cubeAboveFixedPlanCstrs_[i].LB(out);
  }
  else if(iMaxCstrFixPlan_ <= i && i < iMaxCstrS2Plan_)
  {
    auto iPlan = i - cubeAboveFixedPlanCstrs_.size();
    auto iCubeAbove = plans_[iPlan].cubeAbove();
    auto iCubeBelow = plans_[iPlan].cubeBelow();
    cubeAbovePlanFcts_[iCubeAbove].LB(out.head(8));
    cubeAbovePlanFcts_[iCubeBelow].LB(out.tail(8));
  }
  else if(iMaxCstrS2Plan_ <= i && i < iMaxCstrQuatNorm1_)
  {
    out(0) = 1;
  }
  else if(iMaxCstrQuatNorm1_ <= i && i < iMaxCstrVecNorm1_)
  {
    out(0) = 1;
  }
  else
  {
    std::cerr << "WOOOPS! Wrong csrt index " << i << std::endl;
  }
}
void CubeStackProblemOnR::getNonLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  if (i < iMaxCstrFixPlan_)
  {
    cubeAboveFixedPlanCstrs_[i].UB(out);
  }
  else if(iMaxCstrFixPlan_ <= i && i < iMaxCstrS2Plan_)
  {
    auto iPlan = i - cubeAboveFixedPlanCstrs_.size();
    auto iCubeAbove = plans_[iPlan].cubeAbove();
    auto iCubeBelow = plans_[iPlan].cubeBelow();
    cubeAbovePlanFcts_[iCubeAbove].UB(out.head(8));
    cubeAbovePlanFcts_[iCubeBelow].UB(out.tail(8));
  }
  else if(iMaxCstrS2Plan_ <= i && i < iMaxCstrQuatNorm1_)
  {
    out(0) = 1;
  }
  else if(iMaxCstrQuatNorm1_ <= i && i < iMaxCstrVecNorm1_)
  {
    out(0) = 1;
  }
  else
  {
    std::cerr << "WOOOPS! Wrong csrt index " << i << std::endl;
  }
}

size_t CubeStackProblemOnR::numberOfCstr() const
{
  size_t nb = cubeAboveFixedPlanCstrs_.size();
  nb += nPlans_;
  nb += nCubes_;
  nb += nPlans_;
  return nb;
}

Index CubeStackProblemOnR::linCstrDim(size_t) const { return 0; }

Index CubeStackProblemOnR::nonLinCstrDim(size_t i) const
{
  if (i < iMaxCstrFixPlan_)
    return 8;
  else if(iMaxCstrFixPlan_ <= i && i < iMaxCstrS2Plan_)
    return 16;
  else if(iMaxCstrS2Plan_ <= i && i < iMaxCstrQuatNorm1_)
    return 1;
  else if(iMaxCstrQuatNorm1_ <= i && i < iMaxCstrVecNorm1_)
    return 1;
  else
  {
    std::cerr << "WOOOPS! Wrong csrt index " << i << std::endl;
    return 0;
  }
}

std::string CubeStackProblemOnR::getCstrName(const size_t i) const
{
  return cstrNames_[i];
}

void CubeStackProblemOnR::fileForMatlab(std::string fileName,
                                       const Point& x) const
{
  std::cout << "Filing for Matlab" << std::endl;
  std::cout << "x: " << x << std::endl;
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
    std::cout << "Cube " << i << std::endl;
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
