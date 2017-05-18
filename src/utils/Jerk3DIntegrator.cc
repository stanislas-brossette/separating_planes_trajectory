#include <feet-trajectory/utils/Jerk3DIntegrator.hh>
#include <feet-trajectory/utils/MpcCondense.hh>
#include <Eigen/SVD>

namespace feettrajectory
{
void printMat(Eigen::MatrixXd M, std::string Mname)
{
  Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
                               ", ", ", ", "", "", " << ", ";");
  std::cout << Mname << ".resize(" << M.rows() << ", " << M.cols() << ");"
            << std::endl;
  std::cout << Mname << M.format(CommaInitFmt) << std::endl;
}

Jerk3DIntegrator::Jerk3DIntegrator(const double& T, const unsigned long& nIter)
    : T_(T), nIter_(nIter)
{
  Eigen::Matrix3d I3(Eigen::Matrix3d::Identity());
  A_.setZero();
  A_.block(0, 0, 3, 3) = I3;
  A_.block(3, 3, 3, 3) = I3;
  A_.block(6, 6, 3, 3) = I3;
  A_.block(0, 3, 3, 3) = T_ * I3;
  A_.block(3, 6, 3, 3) = T_ * I3;
  A_.block(0, 6, 3, 3) = 0.5 * T_ * T_ * I3;

  B_.setZero();
  B_.topRows(3) = (1. / 6.) * T_ * T_ * T_ * I3;
  B_.middleRows(3, 3) = (1. / 2.) * T_ * T_ * I3;
  B_.bottomRows(3) = T_ * I3;

  mpcCondense(Ux_, Uu_, A_, B_, nIter_);

  SelPos_.resize(3 * nIter_, 9 * nIter_);
  SelVel_.resize(3 * nIter_, 9 * nIter_);
  SelAcc_.resize(3 * nIter_, 9 * nIter_);

  SelPos_.setZero();
  SelVel_.setZero();
  SelAcc_.setZero();

  for (unsigned long i = 0; i < nIter_; i++)
  {
    SelPos_.block(3 * i, 9 * i, 3, 3) = I3;
    SelVel_.block(3 * i, 9 * i + 3, 3, 3) = I3;
    SelAcc_.block(3 * i, 9 * i + 6, 3, 3) = I3;
  }
}

Jerk3DIntegrator::~Jerk3DIntegrator(){};

void Jerk3DIntegrator::getState(Eigen::VectorXd& res, const Eigen::VectorXd& U,
                                const Eigen::VectorXd& X0)
{
  res = Ux_ * X0 + Uu_ * U;
}
void Jerk3DIntegrator::getPos(Eigen::VectorXd& res, const Eigen::VectorXd& U,
                              const Eigen::VectorXd& X0)
{
  res = SelPos_ * (Ux_ * X0 + Uu_ * U);
}
void Jerk3DIntegrator::getVel(Eigen::VectorXd& res, const Eigen::VectorXd& U,
                              const Eigen::VectorXd& X0)
{
  res = SelVel_ * (Ux_ * X0 + Uu_ * U);
}
void Jerk3DIntegrator::getAcc(Eigen::VectorXd& res, const Eigen::VectorXd& U,
                              const Eigen::VectorXd& X0)
{
  res = SelAcc_ * (Ux_ * X0 + Uu_ * U);
}
Eigen::VectorXd Jerk3DIntegrator::getPos(const Eigen::VectorXd& U,
                                const Eigen::VectorXd& X0)
{
  Eigen::VectorXd res(U.size());
  getPos(res, U, X0);
  return res;
}
Eigen::VectorXd Jerk3DIntegrator::getVel(const Eigen::VectorXd& U,
                                const Eigen::VectorXd& X0)
{
  Eigen::VectorXd res(U.size());
  getVel(res, U, X0);
  return res;
}
Eigen::VectorXd Jerk3DIntegrator::getAcc(const Eigen::VectorXd& U,
                                const Eigen::VectorXd& X0)
{
  Eigen::VectorXd res(U.size());
  getAcc(res, U, X0);
  return res;
}
Eigen::VectorXd Jerk3DIntegrator::getState(const Eigen::VectorXd& U,
                                const Eigen::VectorXd& X0)
{
  Eigen::VectorXd res(U.size());
  getState(res, U, X0);
  return res;
}

void Jerk3DIntegrator::jerkFromPos(RefVec jerk, ConstRefVec pos,
                                   ConstRefVec X0, const double& damping) const
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  A = SelPos_*Uu_;
  b = pos - SelPos_ * Ux_ * X0;
  if (damping == 0)
    jerk = A.colPivHouseholderQr().solve(b);
  else
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long nSingularValues = svd.singularValues().size();
    Eigen::VectorXd invS(nSingularValues);
    for (long i = 0; i < nSingularValues; i++)
    {
      invS(i) =
          1 / (svd.singularValues()(i) * svd.singularValues()(i) + damping);
    }
    jerk = svd.matrixV() * Eigen::DiagonalMatrix<double, Eigen::Dynamic>(invS) *
           svd.matrixV().transpose() * A.transpose() * b;
  }
}

} /* feettrajectory */
