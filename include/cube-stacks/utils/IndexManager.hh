#pragma once
#include <iostream>
#include <Eigen/Core>

namespace cubestacks
{
class IndexManager
{
 public:
  IndexManager(const int& nCubes);
  virtual ~IndexManager();
  Eigen::Vector3d getCubeTrans(const int& i, const Eigen::VectorXd& x) const;
  Eigen::Vector4d getCubeQuat(const int& i, const Eigen::VectorXd& x) const;
  double getPlanDist(const int& i, const Eigen::VectorXd& x) const;
  Eigen::Vector3d getPlanNormal(const int& i, const Eigen::VectorXd& x) const;
  int nCubes() const;
  int nPlanes() const;
  int totalDim() const;
  
  int getCubeTransBegin(const int& i) const;
  int getCubeQuatBegin(const int& i) const;
  int getPlaneDistBegin(const int& i) const;
  int getPlaneNormalBegin(const int& i) const;

 private:
  int nCubes_;
  int nPlanes_;
  int totalDim_;
};
} /* cubestacks */
