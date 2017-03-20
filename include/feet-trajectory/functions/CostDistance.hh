#pragma once
#include <iostream>
#include <Eigen/Core>
#include <feet-trajectory/utils/Box.hh>
namespace feettrajectory
{
  class CostDistance
  {
  public:
    CostDistance (const long& nMobileBoxes, const Eigen::Vector3d& initPos, const Eigen::Vector3d& finalPos);
    virtual ~CostDistance ();
    double compute(const Eigen::VectorXd& x);

    const Eigen::MatrixXd& Q() const { return Q_; }
    const Eigen::VectorXd& c() const { return c_; }

  private:
    long n_;
    Eigen::MatrixXd Q_;
    Eigen::VectorXd c_;
  };
} /* feettrajectory */ 
