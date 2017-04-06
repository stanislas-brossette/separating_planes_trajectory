#include <limits>
#include <feet-trajectory/utils/QP.hh>

namespace feettrajectory
{
QP::QP() : QP(0, 0) {}

QP::QP(const Index& dimVar, const Index& dimCstr)
{
  setDimensions(dimVar, dimCstr);
}

void QP::setDimensions(const Index& dimVar, const Index& dimCstr)
{
  dimVar_ = dimVar;
  dimCstr_ = dimCstr;
  A_.resize(dimVar_, dimVar_);
  C_.resize(dimCstr_, dimVar_);
  c_.resize(dimVar_);
  lVar_.resize(dimVar_);
  uVar_.resize(dimVar_);
  l_.resize(dimCstr_);
  u_.resize(dimCstr_);
  A_.setZero();
  C_.setZero();
  c_.setZero();
  l_.setZero();
  u_.setConstant(std::numeric_limits<double>::infinity());
  lVar_.setConstant(-std::numeric_limits<double>::infinity());
  uVar_.setConstant(std::numeric_limits<double>::infinity());
}

QP::~QP() {}

std::ostream& QP::print(std::ostream& o, const Eigen::IOFormat& fmt) const
{
  o << "QP:\n";
  o << "dimVar_:" << dimVar_ << "\n";
  o << "dimCstr_:" << dimCstr_ << "\n";
  o << "A_ quadratic cost:\n" << A_.format(fmt) << "\n";
  o << "c_ linear cost:\n" << c_.transpose().format(fmt) << "\n";
  o << "C_ constraint:\n" << C_.format(fmt) << "\n";
  o << "lVar_ lower bound on variables:\n" << lVar_.transpose().format(fmt) << "\n";
  o << "uVar_ upper bound on variables:\n" << uVar_.transpose().format(fmt) << std::endl;
  o << "l_ lower bound on constraints:\n" << l_.transpose().format(fmt) << "\n";
  o << "u_ upper bound on constraints:\n" << u_.transpose().format(fmt) << std::endl;
  return o;
}

std::ostream& operator<<(std::ostream& o, const QP& f) { return f.print(o); }

} /* feettrajectory */
