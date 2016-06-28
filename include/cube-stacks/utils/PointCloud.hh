#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace cubestacks
{
class PointCloud
{
 public:
  PointCloud(const int& index, const Eigen::VectorXd& v) : index_(index)
  {
    assert(v.size()%3==0 && "Point cloud vector size not a multiple of 3");
    nPoints_ = v.size()/3;
    vertex_.resize(nPoints_);
    for (int i = 0; i < nPoints_; i++) 
    {
      vertex_[i] << v.segment(3*i, 3);
    }
  }

  const std::vector<Eigen::Vector3d>& vertex() const { return vertex_; };
  const Eigen::Vector3d& vertex(long i) const
  {
    return vertex_[static_cast<size_t>(i)];
  }
  const Eigen::Vector3d& vertex(size_t i) const { return vertex_[i]; };
  const int& index() const { return index_;};

 private:
  std::vector<Eigen::Vector3d> vertex_;
  int index_;
  double l_;
};

} /* cubestacks */
