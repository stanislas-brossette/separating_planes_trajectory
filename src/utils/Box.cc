#include <feet-trajectory/utils/Box.hh>

namespace feettrajectory
{
Box::Box(int index, double lx, double ly, double lz, double cx, double cy,
         double cz, bool fixed)
    : index_(index), fixed_(fixed)
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

Box::Box(int index, Eigen::Vector3d size)
    : Box(index, size[0], size[1], size[2])
{
}

Box::Box(int index, Eigen::Vector3d size, Eigen::Vector3d center, bool fix)
    : Box(index, size[0], size[1], size[2], center[0], center[1], center[2],
          fix)
{
}

std::ostream& Box::print(std::ostream& o) const
{
  o << "Box " << index() << std::endl;
  o << "fixed: " << fixed() << std::endl;
  o << "center: " << center().transpose() << std::endl;
  o << "size: " << size().transpose() << std::endl;
  for (size_t i = 0; i < 8; i++)
  {
    o << "vertex(" << i << "): " << vertex(i).transpose() << std::endl;
  }
  return o;
}

std::ostream& operator<<(std::ostream& o, const Box& f) { return f.print(o); }

} /* feettrajectory */
