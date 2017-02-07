#pragma once
#include <iostream>
#include <Eigen/Core>

namespace feettrajectory
{
class IndexManager
{
 public:
  IndexManager(const int& nBoxes);
  virtual ~IndexManager();
  Eigen::Vector3d getBoxTrans(const int& i, const Eigen::VectorXd& x) const;
  Eigen::Vector4d getBoxQuat(const int& i, const Eigen::VectorXd& x) const;
  double getPlanDist(const int& i, const Eigen::VectorXd& x) const;
  Eigen::Vector3d getPlanNormal(const int& i, const Eigen::VectorXd& x) const;
  int nBoxes() const;
  int nPlanes() const;
  int totalDim() const;
  
  int getBoxTransBegin(const int& i) const;
  int getBoxQuatBegin(const int& i) const;
  int getPlaneDistBegin(const int& i) const;
  int getPlaneNormalBegin(const int& i) const;

 private:
  int nBoxes_;
  int nPlanes_;
  int totalDim_;
};
} /* feettrajectory */
