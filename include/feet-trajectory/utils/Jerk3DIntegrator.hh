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
  void getState(Eigen::VectorXd& res, const Eigen::VectorXd& U, const Eigen::VectorXd& X0);
  void getPos(Eigen::VectorXd& res, const Eigen::VectorXd& U, const Eigen::VectorXd& X0);
  void getVel(Eigen::VectorXd& res, const Eigen::VectorXd& U, const Eigen::VectorXd& X0);
  void getAcc(Eigen::VectorXd& res, const Eigen::VectorXd& U, const Eigen::VectorXd& X0);

  const Eigen::MatrixXd& Ux() const { return Ux_; }
  const Eigen::MatrixXd& Uu() const { return Uu_; }
  const Eigen::MatrixXd& SelPos() const { return SelPos_; }
  const Eigen::MatrixXd& SelVel() const { return SelVel_; }
  const Eigen::MatrixXd& SelAcc() const { return SelAcc_; }

  void jerkFromPos(RefVec jerk, ConstRefVec pos, ConstRefVec x0) const;

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
