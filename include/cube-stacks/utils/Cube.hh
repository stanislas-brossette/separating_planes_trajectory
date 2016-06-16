#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace cubestacks
{
class Cube
{
 public:
  Cube(int index, double l = 1.0) : index_(index), l_(l)
  {
    vertex_.resize(8);
    vertex_[0] << -l, -l, -l;
    vertex_[1] << -l, -l, l;
    vertex_[2] << -l, l, -l;
    vertex_[3] << -l, l, l;
    vertex_[4] << l, -l, -l;
    vertex_[5] << l, -l, l;
    vertex_[6] << l, l, -l;
    vertex_[7] << l, l, l;
    vertexMat_.col(0) = vertex_[0];
    vertexMat_.col(1) = vertex_[1];
    vertexMat_.col(2) = vertex_[2];
    vertexMat_.col(3) = vertex_[3];
    vertexMat_.col(4) = vertex_[4];
    vertexMat_.col(5) = vertex_[5];
    vertexMat_.col(6) = vertex_[6];
    vertexMat_.col(7) = vertex_[7];
  }

  const std::vector<Eigen::Vector3d>& vertex() const { return vertex_; };
  const Eigen::Vector3d& vertex(long i) const
  {
    return vertex_[static_cast<size_t>(i)];
  }
  const Eigen::Matrix<double, 3, 8>& vertexMat() const { return vertexMat_; };
  const Eigen::Vector3d& vertex(size_t i) const { return vertex_[i]; };
  const double& l() const { return l_; };
  const int& index() const { return index_;};

 private:
  std::vector<Eigen::Vector3d> vertex_;
  Eigen::Matrix<double, 3, 8> vertexMat_;
  int index_;
  double l_;
};

} /* cubestacks */
