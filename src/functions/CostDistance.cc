#include <feet-trajectory/functions/CostDistance.hh>

namespace feettrajectory
{
CostDistance::CostDistance(const long& nMobileBoxes,
                           ConstRefVec3d initPos,
                           ConstRefVec3d finalPos)
    : n_(3 * nMobileBoxes), Q_(n_, n_), c_(n_)
{
  Q_.setZero();
  c_.setZero();
  c_.head(3) = -2 * initPos;
  c_.tail(3) = -2 * finalPos;
  for (long i = 0; i < n_; ++i)
  {
    for (long j = 0; j < n_; ++j)
    {
      if (j == i)
        Q_(i, j) = 4;
      else if (j == i + 3 || j == i - 3)
        Q_(i, j) = -2;
    }
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
