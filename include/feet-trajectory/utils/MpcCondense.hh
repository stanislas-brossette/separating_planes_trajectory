#pragma once
#include <iostream>
#include <Eigen/Core>

namespace feettrajectory
{
/// @brief Utility class to handle a triple integrator
class MpcCondense
{
 public:
  /// @brief Constructor of the triple integrator
  ///
  /// @param A Matrix multiplying the state (x_k) in the discrete description of
  /// the system
  /// @param B Matrix multiplying the control (u_k) in the discrete description
  /// of the system
  /// @param n number of intervals
  MpcCondense(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
                   const size_t& n);
  virtual ~MpcCondense();
  const Eigen::MatrixXd& Ux() const;
  const Eigen::MatrixXd& Uu() const;

 private:
  /// @brief Dimension of the control variable
  unsigned long dimU_;
  /// @brief Dimension of the state variable
  unsigned long dimX_;
  /// @brief U_x matrix multiplying the x_o term in the condensed form
  Eigen::MatrixXd Ux_;
  /// @brief U_u matrix multiplying the control term in the condensed form
  Eigen::MatrixXd Uu_;
};

} /* feettrajectory */
