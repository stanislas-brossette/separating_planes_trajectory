#include <iostream>
#include <cube-stacks/utils/IndexManager.hh>

namespace cubestacks
{
IndexManager::IndexManager(const int& nCubes) : nCubes_(nCubes)
{
  nPlanes_ = nCubes_ * (nCubes_ - 1) / 2;
  totalDim_ = 7 * nCubes_ + 4 * nPlanes_;
}

IndexManager::~IndexManager() {}

Eigen::Vector3d IndexManager::getCubeTrans(const int& i,
                                           const Eigen::VectorXd& x) const
{
  Eigen::Vector3d res = x.segment(7 * i, 3);
  return res;
}

Eigen::Vector4d IndexManager::getCubeQuat(const int& i,
                                          const Eigen::VectorXd& x) const
{
  Eigen::Vector4d res = x.segment(7 * i + 3, 4);
  return res;
}

double IndexManager::getPlanDist(const int& i, const Eigen::VectorXd& x) const
{
  double res = x(nCubes_ * 7 + i * 4);
  return res;
}

Eigen::Vector3d IndexManager::getPlanNormal(const int& i,
                                            const Eigen::VectorXd& x) const
{
  Eigen::Vector3d res = x.segment(nCubes_ * 7 + i * 4 + 1, 3);
  return res;
}

int IndexManager::nCubes() const { return nCubes_; }

int IndexManager::nPlanes() const { return nPlanes_; }

int IndexManager::totalDim() const { return totalDim_; }

int IndexManager::getCubeTransBegin(const int& i) const { return 7 * i; }

int IndexManager::getCubeQuatBegin(const int& i) const { return 7 * i + 3; }

int IndexManager::getPlaneDistBegin(const int& i) const
{
  return 7 * nCubes_ + 4 * i;
}

int IndexManager::getPlaneNormalBegin(const int& i) const
{
  return 7 * nCubes_ + 4 * i + 1;
}
} /* cubestacks */
