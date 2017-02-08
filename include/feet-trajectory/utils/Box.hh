#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace feettrajectory
{
class Box
{
 public:
  Box(int index, double lx, double ly, double lz, double cx = 0.,
      double cy = 0., double cz = 0.)
      : index_(index)
  {
    size_ << lx, ly, lz;
    center_ << cx, cy, cz;
    vertex_.resize(8);
    vertex_[0] << -lx / 2, -ly / 2, -lz / 2;
    vertex_[1] << -lx / 2, -ly / 2, lz / 2;
    vertex_[2] << -lx / 2, ly / 2, -lz / 2;
    vertex_[3] << -lx / 2, ly / 2, lz / 2;
    vertex_[4] << lx / 2, -ly / 2, -lz / 2;
    vertex_[5] << lx / 2, -ly / 2, lz / 2;
    vertex_[6] << lx / 2, ly / 2, -lz / 2;
    vertex_[7] << lx / 2, ly / 2, lz / 2;
    vertexMat_.col(0) = vertex_[0];
    vertexMat_.col(1) = vertex_[1];
    vertexMat_.col(2) = vertex_[2];
    vertexMat_.col(3) = vertex_[3];
    vertexMat_.col(4) = vertex_[4];
    vertexMat_.col(5) = vertex_[5];
    vertexMat_.col(6) = vertex_[6];
    vertexMat_.col(7) = vertex_[7];
  }

  Box(int index, Eigen::Vector3d size) : Box(index, size[0], size[1], size[2])
  {
  }
  Box(int index, Eigen::Vector3d size, Eigen::Vector3d center)
      : Box(index, size[0], size[1], size[2], center[0], center[1], center[2])
  {
  }

  const std::vector<Eigen::Vector3d>& vertex() const { return vertex_; };
  const Eigen::Matrix<double, 3, 8>& vertexMat() const { return vertexMat_; };
  const Eigen::Vector3d& vertex(size_t i) const { return vertex_[i]; };
  const double& lx() const { return size_[0]; };
  const double& ly() const { return size_[1]; };
  const double& lz() const { return size_[2]; };
  const double& cx() const { return center_[0]; };
  const double& cy() const { return center_[1]; };
  const double& cz() const { return center_[2]; };
  const int& index() const { return index_; };
  const Eigen::Vector3d& center() const { return center_; }
  const Eigen::Vector3d& size() const { return size_; }

  void print() const
  {
    std::cout << "Box " << index() << std::endl;
    std::cout << "center: " << center() << std::endl;
    std::cout << "size: " << size() << std::endl;
    for (size_t i = 0; i < 8; i++)
    {
      std::cout << "vertex(i): " << vertex(i).transpose() << std::endl;
    }
  }

 private:
  std::vector<Eigen::Vector3d> vertex_;
  Eigen::Matrix<double, 3, 8> vertexMat_;
  int index_;
  // Size
  Eigen::Vector3d size_;
  // Center
  Eigen::Vector3d center_;
};

} /* cubestacks */
