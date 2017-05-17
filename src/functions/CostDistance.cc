#include <feet-trajectory/functions/CostDistance.hh>

namespace feettrajectory
{
CostDistance::CostDistance(const long& nMobileBoxes, ConstRefVec3d initPos,
                           ConstRefVec3d finalPos)
    : n_(3 * nMobileBoxes), Q_(n_, n_), c_(n_)
{
  Eigen::Matrix3d I3;
  I3.setIdentity();
  Q_.setZero();
  c_.setZero();
  c_.head(3) = -2 * initPos;
  Q_.block(0, 0, 3, 3) += 2 * I3;
  for (long i = 0; i < nMobileBoxes - 1; ++i)
  {
    Q_.block(3 * i, 3 * i, 3, 3) += 2 * I3;
    Q_.block(3 * i + 3, 3 * i + 3, 3, 3) += 2 * I3;
    Q_.block(3 * i, 3 * i + 3, 3, 3) -= 2 * I3;
    Q_.block(3 * i + 3, 3 * i, 3, 3) -= 2 * I3;
  }
}

CostDistance::~CostDistance() {}

double CostDistance::compute(ConstRefVec x)
{
  assert(x.size() == n_ &&
         "x should be the same size as 3*numberOfMobileBoxes");

  double res = x.dot(Q_ * x) + c_.dot(x);
  return res;
}

void CostDistance::fillQuadCost(RefMat Q, RefVec c) const
{
  Q << Q_;
  c << c_;
}

} /* feettrajectory */
