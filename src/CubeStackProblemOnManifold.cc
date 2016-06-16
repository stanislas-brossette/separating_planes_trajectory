#include <iostream>
#include <fstream>
#include <limits>

#include <Eigen/Geometry>

#include <cube-stacks/CubeStackProblemOnManifold.hh>
#include <cube-stacks/utils/ProblemConfig.hh>

#include <pgsolver/utils/usingManifold.h>
using mnf::RealSpace;
using mnf::SO3;
using mnf::S2;
using mnf::CartesianProduct;
using mnf::ExpMapMatrix;
namespace cubestacks
{
CubeStackProblemOnManifold::CubeStackProblemOnManifold(
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

  //nCubes_ = static_cast<size_t>(M(0).tangentDim()) / 6;
  for (int i = 0; i < static_cast<int>(nCubes_); ++i)
  {
    Cube c(i, 0.5);
    cubes_.push_back(c);
    cubeAboveFixedPlans_.push_back(CubeAboveFixedPlan(c,normalZPlus_,distZPlus_));
  }

  nPlans_ = (nCubes_ * (nCubes_ - 1)) / 2;

  std::cout << "nCubes_ = \n" << nCubes_ << std::endl;
  std::cout << "nPlans_ = \n" << nPlans_ << std::endl;

  //for (size_t i = 0; i < nCubes_; ++i)
  //{
    //for (size_t j = i + 1; j < nCubes_; ++j)
    //{
      //plans_.push_back(Plan(i, j));
    //}
  //}

  std::stringstream sstm;
  sstm << "Stack" << nCubes_ << "CubesS2";
  name() = sstm.str();

  outRepObjDiff_.resize(1, M.representationDim());
  outRepObjDiff_.setZero();
  outRepLinCstrDiff_.resize(static_cast<Index>(8 * nCubes_),
                            M.representationDim());
  outRepLinCstrDiff_.setZero();
  outRepNonLinCstrDiff_.resize(static_cast<Index>(16 * nPlans_),
                               M.representationDim());
  outRepNonLinCstrDiff_.setZero();
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
  SO3<ExpMapMatrix>* SO3_ = new SO3<ExpMapMatrix>();
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

  for (int i = 0; i < out.size(); i++) 
    out(i) = -infinity;
}
void CubeStackProblemOnManifold::getTangentUB(RefVec out) const
{
  assert(out.size() == M().tangentDim() && "wrong size");
  double infinity = std::numeric_limits<double>::infinity();

  for (int i = 0; i < out.size(); i++) 
    out(i) = -infinity;
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
        0, 0, 0, 0, 0, 0;
  }
  M().applyDiffRetractation(out, outRepObjDiff_, x().value());
}

void CubeStackProblemOnManifold::evalLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
  //if (i == 0)
  //{
    //// Constraint 0
    //// All points of all the cubes above the ground

    //out.setZero();
    //// Eigen::Vector3d groundNormal(0, 0, 1);
    //// All cubes above ground
    //cstrAllCubeAboveFixedPlan(out, phi_x_z()(0), -wallZPlus_, groundNormal);
  //}
  //else if (i == 1)
  //{
    //// Constraint 1
    //// This constraint does not have a linear part
  //}
  //else if (i == 2)
  //{
    //// Constraint 2
    //// Front X wall
    //out.setZero();
    //// Eigen::Vector3d normalXPlus(-1, 0, 0);
    //// All cubes above a plan
    //cstrAllCubeAboveFixedPlan(out, phi_x_z()(0), wallXPlus_, normalXPlus);
  //}
  //else if (i == 3)
  //{
    //// Constraint 2
    //// Back X wall
    //out.setZero();
    //// Eigen::Vector3d normalXMinus(1, 0, 0);
    //// All cubes above a plan
    //cstrAllCubeAboveFixedPlan(out, phi_x_z()(0), wallXMinus_, normalXMinus);
  //}
  //else if (i == 4)
  //{
    //// Constraint 4
    //// Front Y wall
    //out.setZero();
    //// Eigen::Vector3d normalYPlus(0, -1, 0);
    //// All cubes above a plan
    //cstrAllCubeAboveFixedPlan(out, phi_x_z()(0), wallYPlus_, normalYPlus);
  //}
  //else if (i == 5)
  //{
    //// Constraint 5
    //// Back Y wall
    //out.setZero();
    //// Eigen::Vector3d normalYMinus(0, 1, 0);
    //// All cubes above a plan
    //cstrAllCubeAboveFixedPlan(out, phi_x_z()(0), wallYMinus_, normalYMinus);
  //}
}
void CubeStackProblemOnManifold::evalLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == linCstrDim(i) && "Wrong rows size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");

  //if (i == 0)
  //{
    //// Constraint 0
    //// All points of the cube above the ground
    //outRepLinCstrDiff_.setZero();
    //Eigen::Vector3d groundNormal(0, 0, 1);
    //cstrAllCubeAboveFixedPlanDiff(outRepLinCstrDiff_, groundNormal);
    //M().applyDiffRetractation(out, outRepLinCstrDiff_, x().value());
  //}
  //else if (i == 1)
  //{
    //// Constraint 1
    //// This constraint does not have a linear part
  //}
  //else if (i == 2)
  //{
    //// Constraint 0
    //// All points of the cube above the ground
    //outRepLinCstrDiff_.setZero();
    //// Eigen::Vector3d normalXPlus(-1, 0, 0);
    //cstrAllCubeAboveFixedPlanDiff(outRepLinCstrDiff_, normalXPlus);
    //M().applyDiffRetractation(out, outRepLinCstrDiff_, x().value());
  //}
  //else if (i == 3)
  //{
    //// Constraint 0
    //// All points of the cube above the ground
    //outRepLinCstrDiff_.setZero();
    //// Eigen::Vector3d normalXMinus(1, 0, 0);
    //cstrAllCubeAboveFixedPlanDiff(outRepLinCstrDiff_, normalXMinus);
    //M().applyDiffRetractation(out, outRepLinCstrDiff_, x().value());
  //}
  //else if (i == 4)
  //{
    //// Constraint 0
    //// All points of the cube above the ground
    //outRepLinCstrDiff_.setZero();
    //// Eigen::Vector3d normalYPlus(0, -1, 0);
    //cstrAllCubeAboveFixedPlanDiff(outRepLinCstrDiff_, normalYPlus);
    //M().applyDiffRetractation(out, outRepLinCstrDiff_, x().value());
  //}
  //else if (i == 5)
  //{
    //// Constraint 0
    //// All points of the cube above the ground
    //outRepLinCstrDiff_.setZero();
    //// Eigen::Vector3d normalYMinus(0, 1, 0);
    //cstrAllCubeAboveFixedPlanDiff(outRepLinCstrDiff_, normalYMinus);
    //M().applyDiffRetractation(out, outRepLinCstrDiff_, x().value());
  //}
}
void CubeStackProblemOnManifold::getLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
  //if (i == 0)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << 0, 0, 0, 0, 0, 0, 0, 0;
    //}
  //}
  //else if (i == 1)
  //{
    //// Constraint 1
    //// This constraint does not have a linear part
  //}
  //else if (i == 2)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << 0, 0, 0, 0, 0, 0, 0, 0;
    //}
  //}
  //else if (i == 3)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << 0, 0, 0, 0, 0, 0, 0, 0;
    //}
  //}
  //else if (i == 4)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << 0, 0, 0, 0, 0, 0, 0, 0;
    //}
  //}
  //else if (i == 5)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << 0, 0, 0, 0, 0, 0, 0, 0;
    //}
  //}
}
void CubeStackProblemOnManifold::getLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == linCstrDim(i) && "wrong size");
  //double inf = std::numeric_limits<double>::infinity();
  //if (i == 0)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << inf, inf, inf, inf, inf, inf, inf, inf;
    //}
  //}
  //else if (i == 1)
  //{
    //// Constraint 1
    //// This constraint does not have a linear part
  //}
  //else if (i == 2)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << inf, inf, inf, inf, inf, inf, inf, inf;
    //}
  //}
  //else if (i == 3)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << inf, inf, inf, inf, inf, inf, inf, inf;
    //}
  //}
  //else if (i == 4)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << inf, inf, inf, inf, inf, inf, inf, inf;
    //}
  //}
  //else if (i == 5)
  //{
    //for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
    //{
      //out.block(8 * i, 0, 8, 1) << inf, inf, inf, inf, inf, inf, inf, inf;
    //}
  //}
}

void CubeStackProblemOnManifold::evalNonLinCstr(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  out.setZero();

  if(i < cubeAboveFixedPlans_.size())
  {
    unsigned long iC = static_cast<unsigned long>(cubeAboveFixedPlans_[i].cube().index());
    Eigen::Vector3d trans = phi_x_z()(0)(iC)[0];
    Eigen::Vector4d quat = phi_x_z()(0)(iC)[1];
    cubeAboveFixedPlans_[i].compute(out.block(8 * static_cast<long>(i), 0, 8, 1), trans, quat); 
  }

  //if (i == 0)
  //{
    //for (size_t i = 0; i < cubes_.size(); i++) 
    //{
      
    //}
    //// Constraint 0
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 1)
  //{
    //// Constraint 1
    //// All planes between cubes
    //out.setZero();
    //Eigen::Vector3d vertex;
    //Eigen::Matrix<double, 8, 1> res;
    //res.setZero();
    //Eigen::Vector3d normalP, posC;
    //for (size_t iP = 0; iP < nPlans_; ++iP)
    //{
      //Eigen::Map<const Eigen::Vector3d> normalP(phi_x_z()(1)(iP)[1].data());
      //double distP = phi_x_z()(1)(iP)[0][0];
      //// Cube Above
      //size_t iC = plans_[iP].cubeAbove();
      //posC = phi_x_z()(0)(iC)[0];
      //Eigen::Map<const Eigen::Matrix3d> rotCA(phi_x_z()(0)(iC)[1].data());
      //for (Index j = 0; j < 8; ++j)
      //{
        //vertex = posC + rotCA * (cubes_[iC].vertex(static_cast<size_t>(j)));
        //res(j) = vertex.dot(normalP) - distP;
      //}
      //out.block(static_cast<Index>(iP * 16), 0, 8, 1) = res;
      //res.setZero();
      //// Cube Below
      //iC = plans_[iP].cubeBelow();
      //posC = phi_x_z()(0)(iC)[0];
      //Eigen::Map<const Eigen::Matrix3d> rotCB(phi_x_z()(0)(iC)[1].data());
      //for (Index j = 0; j < 8; ++j)
      //{
        //vertex = posC + rotCB * (cubes_[iC].vertex(static_cast<size_t>(j)));
        //res(j) = -vertex.dot(normalP) + distP;
      //}
      //out.block(static_cast<Index>(iP * 16 + 8), 0, 8, 1) = res;
      //res.setZero();
    //}
  //}
  //else if (i == 2)
  //{
    //// Constraint 2
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 3)
  //{
    //// Constraint 3
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 4)
  //{
    //// Constraint 4
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 5)
  //{
    //// Constraint 5
    //// This constraint does not have a nonlinear part
  //}
}
void CubeStackProblemOnManifold::evalNonLinCstrDiff(RefMat out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.rows() == nonLinCstrDim(i) && "wrong row size");
  assert(out.cols() == M().tangentDim() && "Wrong cols size");

  if(i < cubeAboveFixedPlans_.size())
  {
    int iC = cubeAboveFixedPlans_[i].cube().index();
    Eigen::Vector3d trans = phi_x_z()(0)(static_cast<size_t>(iC))[0];
    Eigen::Vector4d quat = phi_x_z()(0)(static_cast<size_t>(iC))[1];
    cubeAboveFixedPlans_[i].diffTrans(out.block(8 * static_cast<long>(i), 7* iC, 8, 3), trans, quat); 
    cubeAboveFixedPlans_[i].diffQuat(out.block(8 * static_cast<long>(i), 7* iC + 3, 8, 4), trans, quat); 
  }

  //if (i == 0) {
    //// Constraint 0
    //// This constraint does not have a nonlinear part
  //}

  //else if (i == 1)
  //{
    //// Constraint 1
    //// All planes between cubes
    //outRepNonLinCstrDiff_.setZero();
    //Eigen::Matrix<double, 8, 12> outRepBlockCube;
    //Eigen::Matrix<double, 8, 4> outRepBlockPlane;
    //Eigen::Vector3d nP, v, vP, pC;  // normal of the plane, vertex of the cube
                                    //// and vertex of the cube wrt the plane and
                                    //// position of the cube
    //for (size_t iP = 0; iP < nPlans_; ++iP)
    //{
      //Eigen::Map<const Eigen::Vector3d> nP(phi_x_z()(1)(iP)[1].data());
      //// Cube Above
      //size_t iC = plans_[iP].cubeAbove();
      //pC = phi_x_z()(0)(iC)[0];
      //Eigen::Map<const Eigen::Matrix3d> rCA(phi_x_z()(0)(iC)[1].data());
      //for (Index j = 0; j < 8; ++j)
      //{
        //v = cubes_[i].vertex(static_cast<size_t>(j));
        //outRepBlockCube.row(j) << nP(0), nP(1), nP(2), v(0) * nP(0),
            //v(0) * nP(1), v(0) * nP(2), v(1) * nP(0), v(1) * nP(1),
            //v(1) * nP(2), v(2) * nP(0), v(2) * nP(1), v(2) * nP(2);
        //vP = pC + rCA * v;
        //outRepBlockPlane.row(j) << -1, vP(0), vP(1), vP(2);
      //}
      //outRepNonLinCstrDiff_.block(static_cast<Index>(iP) * 16,
                                  //static_cast<Index>(12 * iC), 8, 12) =
          //outRepBlockCube;
      //outRepNonLinCstrDiff_.block(static_cast<Index>(iP) * 16,
                                  //static_cast<Index>(12 * nCubes_ + 4 * iP), 8,
                                  //4) = outRepBlockPlane;
      //// Cube Below
      //iC = plans_[iP].cubeBelow();
      //pC = phi_x_z()(0)(iC)[0];
      //Eigen::Map<const Eigen::Matrix3d> rCB(phi_x_z()(0)(iC)[1].data());
      //for (Index j = 0; j < 8; ++j)
      //{
        //v = cubes_[i].vertex(static_cast<size_t>(j));
        //outRepBlockCube.row(j) << -nP(0), -nP(1), -nP(2), -v(0) * nP(0),
            //-v(0) * nP(1), -v(0) * nP(2), -v(1) * nP(0), -v(1) * nP(1),
            //-v(1) * nP(2), -v(2) * nP(0), -v(2) * nP(1), -v(2) * nP(2);
        //vP = pC + rCB * v;
        //outRepBlockPlane.row(j) << 1, -vP(0), -vP(1), -vP(2);
      //}
      //outRepNonLinCstrDiff_.block(static_cast<Index>(iP) * 16 + 8,
                                  //static_cast<Index>(12 * iC), 8, 12) =
          //outRepBlockCube;
      //outRepNonLinCstrDiff_.block(static_cast<Index>(iP) * 16 + 8,
                                  //static_cast<Index>(12 * nCubes_ + 4 * iP), 8,
                                  //4) = outRepBlockPlane;
    //}
    //M().applyDiffRetractation(out, outRepNonLinCstrDiff_, x().value());
  //}
  //else if (i == 2)
  //{
    //// Constraint 2
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 3)
  //{
    //// Constraint 3
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 4)
  //{
    //// Constraint 4
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 5)
  //{
    //// Constraint 5
    //// This constraint does not have a nonlinear part
  //}
}
void CubeStackProblemOnManifold::getNonLinCstrLB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  if(i < cubeAboveFixedPlans_.size())
  {
    cubeAboveFixedPlans_[i].LB(out.block(8 * static_cast<long>(i), 0, 8, 1)); 
  }
  //if (i == 0)
  //{
    //// Constraint 0
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 1)
  //{
    //// Constraint 1
    //for (Index i = 0; i < static_cast<Index>(nPlans_); ++i)
    //{
      //out.block(16 * i, 0, 16, 1) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          //0, 0;
    //}
  //}
  //else if (i == 2)
  //{
    //// Constraint 2
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 3)
  //{
    //// Constraint 3
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 4)
  //{
    //// Constraint 4
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 5)
  //{
    //// Constraint 5
    //// This constraint does not have a nonlinear part
  //}
}
void CubeStackProblemOnManifold::getNonLinCstrUB(RefVec out, size_t i) const
{
  assert(i < numberOfCstr() && "This constraint index is out of bounds");
  assert(out.size() == nonLinCstrDim(i) && "wrong size");
  //double inf = std::numeric_limits<double>::infinity();
  if(i < cubeAboveFixedPlans_.size())
  {
    cubeAboveFixedPlans_[i].LB(out.block(8 * static_cast<long>(i), 0, 8, 1)); 
  }
  //else if (i == 1)
  //{
    //// Constraint 0
    //for (Index i = 0; i < static_cast<Index>(nPlans_); ++i)
    //{
      //out.block(16 * i, 0, 16, 1) << inf, inf, inf, inf, inf, inf, inf, inf,
          //inf, inf, inf, inf, inf, inf, inf, inf;
    //}
  //}
  //else if (i == 2)
  //{
    //// Constraint 2
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 3)
  //{
    //// Constraint 3
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 4)
  //{
    //// Constraint 4
    //// This constraint does not have a nonlinear part
  //}
  //else if (i == 5)
  //{
    //// Constraint 5
    //// This constraint does not have a nonlinear part
  //}
}

size_t CubeStackProblemOnManifold::numberOfCstr() const
{
  size_t nb = nCubes_;
  return nb;
}

Index CubeStackProblemOnManifold::linCstrDim(size_t i) const
{
  //Index nCubesI = static_cast<Index>(nCubes_);
  Index linDim[] = {0};
  return linDim[static_cast<size_t>(i)];
}

Index CubeStackProblemOnManifold::nonLinCstrDim(size_t i) const
{
  //Index nPlanesI = static_cast<Index>(nPlans_);
  //Index nCubesI = static_cast<Index>(nCubes_);
  Index nonLinDim[] = {8};
  return nonLinDim[static_cast<size_t>(i)];
}

void CubeStackProblemOnManifold::printState() const
{
  // Eigen::Vector3d pos = x()[0];
  // Eigen::Map<const Eigen::Matrix3d> rot(x()[1].data());
  // Eigen::Vector3d vertex;
  // Eigen::Vector3d groundNormal(0, 0, 1);
  // Eigen::VectorXd res(8);
  // std::cout << "Cube:" << std::endl;
  // std::cout << "Center: " << pos.transpose() << std::endl;
  // for(Index i = 0; i < 8; ++i)
  //{
  //  vertex = pos+rot*cubes_[0].vertex()[(size_t)i];
  //  std::cout << "vertex " << i << ": " << vertex.transpose() << std::endl;
  //}
}

void CubeStackProblemOnManifold::cstrAllCubeAboveFixedPlan(
    RefMat out, const ConstSubPoint& pointAllCubes, const double& distPlan,
    const Eigen::Vector3d& normal) const
{
  for (size_t i = 0; i < nCubes_; ++i)
  {
    cstrCubeAboveFixedPlan(out.block<8, 1>(static_cast<Index>(i) * 8, 0),
                           pointAllCubes(i), cubes_[i], distPlan, normal);
  }
}

void CubeStackProblemOnManifold::cstrCubeAboveFixedPlan(
    RefMat out, const ConstSubPoint& pointCube, const Cube& cube,
    const double& distPlan, const Eigen::Vector3d& normal) const
{
  Eigen::Vector3d vertex;
  Eigen::Vector3d posC = pointCube[0];
  Eigen::Map<const Eigen::Matrix3d> rotC(pointCube[1].data());
  for (Index j = 0; j < 8; ++j)
  {
    vertex = posC + rotC * (cube.vertex(static_cast<size_t>(j)));
    out(j, 0) = distPlan + vertex.dot(normal);
  }
}

void CubeStackProblemOnManifold::cstrAllCubeAboveFixedPlanDiff(
    RefMat out, const Eigen::Vector3d& normal) const
{
  for (Index i = 0; i < static_cast<Index>(nCubes_); ++i)
  {
    cstrCubeAboveFixedPlanDiff(out.block(8 * i, 12 * i, 8, 12),
                               cubes_[static_cast<size_t>(i)], normal);
  }
}

void CubeStackProblemOnManifold::cstrCubeAboveFixedPlanDiff(
    RefMat out, const Cube& cube, const Eigen::Vector3d& normal) const
{
  for (Index j = 0; j < 8; ++j)
  {
    Eigen::Vector3d v = cube.vertex(static_cast<size_t>(j));
    out.row(j) << normal(0), normal(1), normal(2), v(0) * normal(0),
        v(0) * normal(1), v(0) * normal(2), v(1) * normal(0), v(1) * normal(1),
        v(1) * normal(2), v(2) * normal(0), v(2) * normal(1), v(2) * normal(2);
  }
}

//void CubeStackProblemOnManifold::fileForMatlab(std::string fileName,
                                       //const Point& x) const
//{
  //std::ofstream myFile;
  //myFile.open(fileName);

  //myFile << "clear all;" << std::endl;
  //myFile << "material shiny;" << std::endl;
  //myFile << "alpha(\'color\');" << std::endl;
  //myFile << "alphamap(\'rampdown\');" << std::endl;
  //myFile << "xlabel(\'x\');" << std::endl;
  //myFile << "ylabel(\'y\');" << std::endl;
  //myFile << "zlabel(\'z\');" << std::endl;
  //myFile << "grid \'on\';" << std::endl;
  //myFile << "view(30,30);" << std::endl;
  //myFile << "axis(\'equal\');   " << std::endl;

  //double minX(-2.0);
  //double maxX(2.0);
  //double minY(-2.0);
  //double maxY(2.0);
  //double minZ(-2.0);
  //double maxZ(2.0);
  //for (size_t i = 0; i < nCubes_; ++i)
  //{
    //Eigen::Vector3d pC = x(0)(i)[0];
    //minX = fmin(minX, pC(0) - cubes_[i].l());
    //maxX = fmax(maxX, pC(0) + cubes_[i].l());
    //minY = fmin(minY, pC(1) - cubes_[i].l());
    //maxY = fmax(maxY, pC(1) + cubes_[i].l());
    //minZ = fmin(minZ, pC(2) - cubes_[i].l());
    //maxZ = fmax(maxZ, pC(2) + cubes_[i].l());
    //Eigen::Map<const Eigen::Matrix3d> rotC(x(0)(i)[1].data());
    //myFile << "vert" << i << " = [ ";
    //for (Index j = 0; j < 8; ++j)
    //{
      //Eigen::Vector3d v =
          //pC + rotC * (cubes_[i].vertex(static_cast<size_t>(j)));
      //myFile << v.transpose();

      //if (j < 7)
      //{
        //myFile << "; ..." << std::endl;
      //}
      //else
      //{
        //myFile << "];" << std::endl;
      //}
    //}
    //myFile << "fac" << i
           //<< " = [1 2 4 3; 2 6 8 4; 4 3 7 8; 1 5 7 3; 1 2 6 5; 5 6 8 7];"
           //<< std::endl;
    //myFile << "patch( \'Faces\', fac" << i << ", \'Vertices\', vert" << i
           //<< ", \'FaceColor\', \'r\');" << std::endl;
  //}
  //for (size_t i = 0; i < nPlans_; ++i)
  //{
    //double r = x(1)(i)[0][0];
    //Eigen::Vector3d n = x(1)(i)[1];
    //size_t cubeBelow = plans_[i].cubeBelow();
    //size_t cubeAbove = plans_[i].cubeAbove();
    //Eigen::Vector3d Cb = x(0)(cubeBelow)[0];
    //Eigen::Vector3d Cu = x(0)(cubeAbove)[0];

    //double lambda = (r - n.dot(Cb)) / n.dot(Cu - Cb);
    //Eigen::Vector3d zeroPlane = Cb + lambda * (Cu - Cb);

    //Eigen::Vector3d t(n(1), -n(0), 0);
    //t = t / t.norm();

    //Eigen::Vector3d b = n.cross(t);
    //b = b / b.norm();

    //myFile << "vert" << i << " = [ ";

    //Eigen::Vector3d vertex;
    //vertex = zeroPlane + t + b;
    //myFile << vertex.transpose() << "; ..." << std::endl;
    //vertex = zeroPlane + t - b;
    //myFile << vertex.transpose() << "; ..." << std::endl;
    //vertex = zeroPlane - t - b;
    //myFile << vertex.transpose() << "; ..." << std::endl;
    //vertex = zeroPlane - t + b;
    //myFile << vertex.transpose() << "; ..." << std::endl;
    //myFile << "];" << std::endl;
    //myFile << "fac" << i << " = [1 2 3 4];" << std::endl;
    //myFile << "patch( \'Faces\', fac" << i << ", \'Vertices\', vert" << i
           //<< ", \'FaceColor\', \'y\', 'FaceAlpha', 0.5);" << std::endl;
  //}
  //myFile << "xlim([" << minX << " " << maxX << "]);" << std::endl;
  //myFile << "ylim([" << minY << " " << maxY << "]);" << std::endl;
  //myFile << "zlim([" << minZ << " " << maxZ << "]);" << std::endl;

  //myFile << "vertPlanZPlus = [ ";
  //myFile << -wallXMinus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << -wallXMinus_ << ", " << wallYPlus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << wallYPlus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ << "];"
         //<< std::endl;
  //myFile << "facZPlus = [1 2 3 4];" << std::endl;
  //myFile << "patch( \'Faces\', facZPlus, \'Vertices\', vertPlanZPlus, "
            //"\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  //myFile << "vertPlanXPlus = [ ";
  //myFile << wallXPlus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << wallYPlus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << wallYPlus_ << ", " << wallZPlus_ + 3 << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ + 3 << "];"
         //<< std::endl;
  //myFile << "facXPlus = [1 2 3 4];" << std::endl;
  //myFile << "patch( \'Faces\', facXPlus, \'Vertices\', vertPlanXPlus, "
            //"\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  //myFile << "vertPlanXMinus = [ ";
  //myFile << -wallXMinus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << -wallXMinus_ << ", " << wallYPlus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << -wallXMinus_ << ", " << wallYPlus_ << ", " << wallZPlus_ + 3
         //<< ";..." << std::endl;
  //myFile << -wallXMinus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ + 3
         //<< "];" << std::endl;
  //myFile << "facXMinus = [1 2 3 4];" << std::endl;
  //myFile << "patch( \'Faces\', facXMinus, \'Vertices\', vertPlanXMinus, "
            //"\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  //myFile << "vertPlanYPlus = [ ";
  //myFile << -wallXMinus_ << ", " << wallYPlus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << wallYPlus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << wallYPlus_ << ", " << wallZPlus_ + 3 << ";..."
         //<< std::endl;
  //myFile << -wallXMinus_ << ", " << wallYPlus_ << ", " << wallZPlus_ + 3 << "];"
         //<< std::endl;
  //myFile << "facYPlus = [1 2 3 4];" << std::endl;
  //myFile << "patch( \'Faces\', facYPlus, \'Vertices\', vertPlanYPlus, "
            //"\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  //myFile << "vertPlanYMinus = [ ";
  //myFile << -wallXMinus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ << ";..."
         //<< std::endl;
  //myFile << wallXPlus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ + 3
         //<< ";..." << std::endl;
  //myFile << -wallXMinus_ << ", " << -wallYMinus_ << ", " << wallZPlus_ + 3
         //<< "];" << std::endl;
  //myFile << "facYMinus = [1 2 3 4];" << std::endl;
  //myFile << "patch( \'Faces\', facYMinus, \'Vertices\', vertPlanYMinus, "
            //"\'FaceColor\', \'b\', 'FaceAlpha', 0.5);" << std::endl;
  //myFile.close();
//}
}
