#include <iostream>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <feet-trajectory/utils/Jerk3DIntegrator.hh>

using namespace feettrajectory;

BOOST_AUTO_TEST_CASE(testJerk3DIntegrator)
{
  Jerk3DIntegrator integ(0.1, 4);
  BOOST_CHECK(true);
}
