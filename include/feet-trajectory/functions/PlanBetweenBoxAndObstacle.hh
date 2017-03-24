#pragma once
#include <iostream>
#include <Eigen/Core>
#include <feet-trajectory/utils/Box.hh>
#include <feet-trajectory/utils/defs.hh>

namespace feettrajectory
{
  class PlanBetweenBoxAndObstacle
  {
  public:
    PlanBetweenBoxAndObstacle ();
    virtual ~PlanBetweenBoxAndObstacle ();
    static void fillLinCstr(const Box& box, const Box& obstacle,
                            const double& planD, const Eigen::Vector3d& planN,
                            RefVec lb, RefMat C);

  private:
    /* data */
  };
} /* feettrajectory */ 
