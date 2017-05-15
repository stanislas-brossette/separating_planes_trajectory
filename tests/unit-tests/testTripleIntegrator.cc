#include <iostream>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <feet-trajectory/utils/TripleIntegrator.hh>

using namespace feettrajectory;
// int main(void)
//{
// std::cout << "HelloWorld" << std::endl;
// return 0;
//}

BOOST_AUTO_TEST_CASE(testTripleIntegrator)
{
  Eigen::Matrix3d A;
  Eigen::Vector3d B;
  double T(3);
  A << 1, T, T* T / 2, 0, 1, T, 0, 0, 1;
  B << T* T* T / 6., T* T / 2., T;
  TripleIntegrator tig(A, B, 3);

  Eigen::Matrix<double, 9, 3> expUx;
  Eigen::Matrix<double, 9, 3> expUu;

  expUx << 1. ,  3. ,  4.5,
           0. ,  1. ,  3.,
           0. ,  0. ,  1.,
           1. ,  6. ,  18,
           0. ,  1. ,  6.,
           0. ,  0. ,  1.,
           1. ,  9. ,  40.5,
           0. ,  1. ,  9.,
           0. ,  0. ,  1.;
  expUu << 4.5, 0, 0,
           4.5, 0, 0,
           3, 0, 0,
           31.5, 4.5, 0,
           13.5, 4.5, 0,
           3, 3, 0,
           85.5,31.5, 4.5,
           22.5,13.5, 4.5,
           3, 3, 3;

  BOOST_CHECK(tig.Ux().isApprox(expUx));
  BOOST_CHECK(tig.Uu().isApprox(expUu));
}
