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
      double cy = 0., double cz = 0., bool fixed = false, bool isVirtual = false);
  Box(int index, Eigen::Vector3d size);
  Box(int index, Eigen::Vector3d size, Eigen::Vector3d center, bool fix = false, bool isVirtual = false);

  const std::vector<Eigen::Vector3d>& vertex() const { return vertex_; };
  const Eigen::Matrix<double, 3, 8>& vertexMat() const { return vertexMat_; };
  const Eigen::Vector3d& vertex(size_t i) const;
  const double& lx() const { return size_[0]; };
  const double& ly() const { return size_[1]; };
  const double& lz() const { return size_[2]; };
  const double& cx() const { return center_[0]; };
  const double& cy() const { return center_[1]; };
  const double& cz() const { return center_[2]; };
  const int& index() const { return index_; };
  const bool& fixed() const { return fixed_; };
  const bool& isVirtual() const { return isVirtual_; };
  void setFixed(bool f) { fixed_ = f; };
  void setVirtual(bool v) { isVirtual_ = v; };
  const Eigen::Vector3d& center() const { return center_; }
  const Eigen::Vector3d& size() const { return size_; }

  /// \brief Print method.
  /// \param o output stream.
  /// \return output stream.
  std::ostream& print (std::ostream& o) const;

 private:
  std::vector<Eigen::Vector3d> vertex_;
  Eigen::Matrix<double, 3, 8> vertexMat_;
  int index_;
  // Size
  Eigen::Vector3d size_;
  // Center
  Eigen::Vector3d center_;
  // Is fixed
  bool fixed_;
  // Is virtual (used only for obstacles
  bool isVirtual_;
};

/// \brief Output stream operator for frames.
/// \param o output stream.
/// \param f Box.
/// \return output stream.
std::ostream& operator<<(std::ostream& o, const Box& f);

} /* feettrajectory */
