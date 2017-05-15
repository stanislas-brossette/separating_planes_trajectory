#include <feet-trajectory/utils/TripleIntegrator.hh>

namespace feettrajectory
{
TripleIntegrator::TripleIntegrator(const Eigen::MatrixXd& A,
                                   const Eigen::MatrixXd& B, const size_t& n)
    : Ux_(n * A.rows(), A.cols()), Uu_(n * B.rows(), n * B.cols())
{
  Ux_.setZero();
  Uu_.setZero();
  Eigen::MatrixXd M(A);
  for (int i = 0; i < n; i++)
  {
    Ux_.middleRows(i * A.rows(), A.rows()) = M;
    M = M * A;
  }

  M = B;
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n - i; j++)
    {
      Uu_.block((i + j) * B.rows(), j * B.cols(), B.rows(), B.cols()) = M;
    }
    M = A * M;
  }
}

TripleIntegrator::~TripleIntegrator() {}

const Eigen::MatrixXd& TripleIntegrator::Ux() const { return Ux_; }
const Eigen::MatrixXd& TripleIntegrator::Uu() const { return Uu_; }

} /* feettrajectory */
