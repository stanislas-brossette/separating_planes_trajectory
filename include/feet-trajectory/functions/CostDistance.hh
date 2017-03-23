#pragma once
#include <iostream>
#include <Eigen/Core>
#include <feet-trajectory/utils/Box.hh>
#include <feet-trajectory/utils/defs.hh>
namespace feettrajectory
{
  class CostDistance
  {
  public:
    CostDistance (const long& nMobileBoxes, ConstRefVec3d initPos, ConstRefVec3d finalPos);
    virtual ~CostDistance ();
    double compute(ConstRefVec x);

    const Eigen::MatrixXd& Q() const { return Q_; }
    const Eigen::VectorXd& c() const { return c_; }

    void fillQuadCost(RefMat Q, RefVec c) const;

  private:
    long n_;
    Eigen::MatrixXd Q_;
    Eigen::VectorXd c_;
  };
} /* feettrajectory */ 
