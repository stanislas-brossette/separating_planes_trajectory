#pragma once

#include <iostream>
#include <feet-trajectory/utils/QP.hh>
#include <feet-trajectory/TrajectoryProblem.hh>

namespace feettrajectory
{
class QPWithNFixed : public QP
{
 public:
  QPWithNFixed(const TrajectoryProblem& pb);
  virtual ~QPWithNFixed();

 private:
};

} /* feettrajectory */

