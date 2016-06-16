CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME cube-stacks)
SET(PROJECT_DESCRIPTION "Compare resolution of cube stack problems with different optimization approaches")
SET(PROJECT_URL "git@gite.lirmm.fr:brossette/cube-stacks.git")
SET(PROJECT_DEBUG_POSTFIX "_d")

SETUP_PROJECT()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")

#########################
# External dependencies
SEARCH_FOR_EIGEN()
SET(BOOST_REQUIRED 1.48)
SET(BOOST_COMPONENTS system unit_test_framework)

ADD_REQUIRED_DEPENDENCY("yaml-cpp >= 0.5")
ADD_REQUIRED_DEPENDENCY("roboptim-core > 3.0")
ADD_REQUIRED_DEPENDENCY(roboptim-core-manifold roboptim-core-manifold_d)
ADD_REQUIRED_DEPENDENCY(manifolds manifolds_d)
ADD_REQUIRED_DEPENDENCY(PGSolver PGSolver_d)

########################

ADD_DEFINITIONS(-DROBOPTIM_DO_NOT_CHECK_ALLOCATION)
ADD_DEFINITIONS(-DLOGS_DATA_DIR="${PROJECT_SOURCE_DIR}/logs")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

# Add dependency towards the library in the pkg-config file.
PKG_CONFIG_APPEND_LIBS(cube-stacks)

SETUP_PROJECT_FINALIZE()
#SETUP_PROJECT_CPACK()
