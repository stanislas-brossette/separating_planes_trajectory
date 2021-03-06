#pragma once

#include <iostream>
#include <Eigen/Core>
#include <feet-trajectory/utils/defs.hh>

namespace feettrajectory
{
class QP
{
 public:
  QP();
  QP(const Index& dimVar, const Index& dimCstr);
  void setDimensions(const Index& dimVar, const Index& dimCstr);
  virtual ~QP();

  //setZero on all the terms of the QP
  void reset();

  const Index& dimVar() const { return dimVar_; }
  const Index& dimCstr() const { return dimCstr_; }
  const Eigen::MatrixXd& A() const { return A_; }
  Eigen::MatrixXd Acopy() const { return A_; }
  const Eigen::MatrixXd& C() const { return C_; }
  const Eigen::VectorXd& c() const { return c_; }
  const Eigen::VectorXd& lVar() const { return lVar_; }
  const Eigen::VectorXd& uVar() const { return uVar_; }
  const Eigen::VectorXd& l() const { return l_; }
  const Eigen::VectorXd& u() const { return u_; }

  /// \brief Print method.
  /// \param o output stream.
  /// \return output stream.
  std::ostream& print (std::ostream& o, const Eigen::IOFormat& fmt = fmt::custom) const;

 protected:
  Index dimVar_;           // Size variables
  Index dimCstr_;          // Size constraints
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
