#include <feet-trajectory/utils/QP.hh>

namespace feettrajectory
{
QP::QP() : QP(0, 0) {}

QP::QP(const long& n, const long& m) { setDimensions(n, m); }

void QP::setDimensions(const long& n, const long& m)
{
  n_ = n;
  m_ = m;
  A_.resize(n_, n_);
  C_.resize(m_, n_);
  c_.resize(n_);
  l_.resize(n_ + m_);
  u_.resize(n_ + m_);
  A_.setZero();
  C_.setZero();
  c_.setZero();
  l_.setZero();
  u_.setZero();
}

QP::~QP() {}
} /* feettrajectory */
