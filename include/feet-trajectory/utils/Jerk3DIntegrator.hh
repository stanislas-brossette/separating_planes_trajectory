#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <feet-trajectory/utils/defs.hh>

namespace feettrajectory
{
/// @brief Class responsible for integrating a 3D jerk control into positions,
/// velocities and accelerations
class Jerk3DIntegrator
{
 public:
  Jerk3DIntegrator(const double& T, const unsigned long& nIter);
  virtual ~Jerk3DIntegrator();
  void getState(RefVec res, ConstRefVec U, ConstRefVec X0) const;
  void getPos(RefVec res, ConstRefVec U, ConstRefVec X0) const;
  void getVel(RefVec res, ConstRefVec U, ConstRefVec X0) const;
  void getAcc(RefVec res, ConstRefVec U, ConstRefVec X0) const;

  Eigen::VectorXd getState(ConstRefVec U, ConstRefVec X0) const;
  Eigen::VectorXd getPos(ConstRefVec U, ConstRefVec X0) const;
  Eigen::VectorXd getVel(ConstRefVec U, ConstRefVec X0) const;
  Eigen::VectorXd getAcc(ConstRefVec U, ConstRefVec X0) const;

  const Eigen::MatrixXd& Ux() const { return Ux_; }
  const Eigen::MatrixXd& Uu() const { return Uu_; }
  const Eigen::MatrixXd& SelPos() const { return SelPos_; }
  const Eigen::MatrixXd& SelVel() const { return SelVel_; }
  const Eigen::MatrixXd& SelAcc() const { return SelAcc_; }

  /// @brief Computes a jerk that minimizes: || Ux.X0 + Uu.jerk - pos ||^2 + damping*jerk^2
  ///
  /// @param jerk result jerk
  /// @param pos desired position
  /// @param x0 initial state
  /// @param damping damping factor
  void jerkFromPos(RefVec jerk, ConstRefVec pos, ConstRefVec x0,
                   const double& damping = 0) const;

 private:
  double T_;  // time step between iterates
  unsigned long nIter_;  // number of iterates
  Eigen::Matrix<double, 9, 9> A_;  // Matrix A such that X_{k+1} = A X_k + B U_k
                                   // with X_k the state and U_k the control
                                   // (here the jerk3D)
  Eigen::Matrix<double, 9, 3> B_;  // Matrix B such that X_{k+1} = A X_k + B U_k
                                   // with X_k the state and U_k the control
                                   // (here the jerk3D)
  Eigen::MatrixXd Ux_;  // Condensed matrix Uu such that Vx = Ux.X0 + Uu.Vu
  Eigen::MatrixXd Uu_;  // Condensed matrix Uu such that Vx = Ux.X0 + Uu.Vu
  Eigen::MatrixXd SelPos_, SelVel_, SelAcc_;  // Selection matrices for
                                              // positions, velocities and
                                              // accelerations respectively
};
} /* feettrajectory */
