#include <iostream>
#include <feet-trajectory/utils/IndexManager.hh>

namespace feettrajectory
{
IndexManager::IndexManager(const int& nBoxes) : nBoxes_(nBoxes)
{
  nPlanes_ = nBoxes_ * (nBoxes_ - 1) / 2;
  totalDim_ = 7 * nBoxes_ + 4 * nPlanes_;
}

IndexManager::~IndexManager() {}

Eigen::Vector3d IndexManager::getBoxTrans(const int& i,
                                          const Eigen::VectorXd& x) const
{
  Eigen::Vector3d res = x.segment(7 * i, 3);
  return res;
}

Eigen::Vector4d IndexManager::getBoxQuat(const int& i,
                                         const Eigen::VectorXd& x) const
{
  Eigen::Vector4d res = x.segment(7 * i + 3, 4);
  return res;
}

double IndexManager::getPlanDist(const int& i, const Eigen::VectorXd& x) const
{
  double res = x(nBoxes_ * 7 + i * 4);
  return res;
}

Eigen::Vector3d IndexManager::getPlanNormal(const int& i,
                                            const Eigen::VectorXd& x) const
{
  Eigen::Vector3d res = x.segment(nBoxes_ * 7 + i * 4 + 1, 3);
  return res;
}

int IndexManager::nBoxes() const { return nBoxes_; }

int IndexManager::nPlanes() const { return nPlanes_; }

int IndexManager::totalDim() const { return totalDim_; }

int IndexManager::getBoxTransBegin(const int& i) const { return 7 * i; }

int IndexManager::getBoxQuatBegin(const int& i) const { return 7 * i + 3; }

int IndexManager::getPlaneDistBegin(const int& i) const
{
  return 7 * nBoxes_ + 4 * i;
}

int IndexManager::getPlaneNormalBegin(const int& i) const
{
  return 7 * nBoxes_ + 4 * i + 1;
}
} /* feettrajectory */
