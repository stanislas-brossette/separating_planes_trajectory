
#include <feet-trajectory/functions/PlanBetweenBoxAndObstacle.hh>
namespace feettrajectory
{
void PlanBetweenBoxAndObstacle::fillLinCstr(const Box& box, const Box& obstacle,
                                            const double& planD,
                                            const Eigen::Vector3d& planN,
                                            RefVec lb, RefMat C,
                                            double securityDistance)
{
  C << planN.transpose();
  const double inf(std::numeric_limits<double>::infinity());
  double maxObs = -inf;
  double minBox = inf;
  double val = inf;

  for (size_t i = 0; i < box.vertex().size(); i++)
  {
    val = box.vertex(i).dot(planN);
    if (val < minBox) minBox = val;
  }

  for (size_t i = 0; i < obstacle.vertex().size(); i++)
  {
    val = obstacle.vertex(i).dot(planN);
    if (val > maxObs) maxObs = val;
  }

  if (!box.fixed())
    lb << maxObs + planN.dot(obstacle.center()) - minBox + 2 * securityDistance;
  else
    lb << maxObs + planN.dot(obstacle.center()) - minBox -
              planN.dot(box.center()) + 2 * securityDistance;
}
} /* feettrajectory */
