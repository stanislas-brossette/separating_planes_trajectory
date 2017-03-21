#pragma once

#include <iostream>
#include <Eigen/Core>

namespace feettrajectory
{
class QP
{
 public:
  QP();
  QP(const long& dimVar, const long& dimCstr);
  void setDimensions(const long& dimVar, const long& dimCstr);
  virtual ~QP();

  /// \brief Print method.
  /// \param o output stream.
  /// \return output stream.
  std::ostream& print (std::ostream& o) const;

 protected:
  long dimVar_;           // Size variables
  long dimCstr_;          // Size constraints
  Eigen::MatrixXd A_;     // Quadratic cost
  Eigen::MatrixXd C_;     // Constraint
  Eigen::VectorXd c_;     // linear cost
  Eigen::VectorXd lVar_;  // lower bound on variables
  Eigen::VectorXd uVar_;  // upper bound on variables
  Eigen::VectorXd l_;     // lower bound on linear constraints
  Eigen::VectorXd u_;     // upper bound on linear constraints
};

/// \brief Output stream operator for frames.
/// \param o output stream.
/// \param f QP.
/// \return output stream.
std::ostream& operator<<(std::ostream& o, const QP& f);
} /* feettrajectory */
