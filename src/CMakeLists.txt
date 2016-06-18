INCLUDE_DIRECTORIES(BEFORE ${Boost_INCLUDE_DIR})

SET(SOURCES
  #CubeStackProblem.cc
  #CubeStackProblem1Cube.cc
  CubeStackProblemOnManifold.cc
  functions/Cost.cc
  functions/CubeAboveFixedPlan.cc
  functions/CubeAbovePlan.cc
  #functions/PlanBetweenCubes.cc
  utils/ProblemConfig.cc
  utils/quat2mat.cc
  )
SET(HEADERS
  #../include/cube-stacks/CubeStackProblem.hh
  #../include/cube-stacks/CubeStackProblem1Cube.hh
  ../include/cube-stacks/CubeStackProblemOnManifold.hh
  ../include/cube-stacks/functions/Cost.hh
  ../include/cube-stacks/functions/CubeAboveFixedPlan.hh
  ../include/cube-stacks/functions/CubeAbovePlan.hh
  #../include/cube-stacks/functions/PlanBetweenCubes.hh
  ../include/cube-stacks/utils/Cube.hh
  ../include/cube-stacks/utils/Plan.hh
  ../include/cube-stacks/utils/ProblemConfig.hh
  ../include/cube-stacks/utils/quat2mat.hh
  )

ADD_LIBRARY(cube-stacks SHARED ${SOURCES} ${HEADERS})

PKG_CONFIG_USE_DEPENDENCY(cube-stacks roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(cube-stacks roboptim-core-manifold)
PKG_CONFIG_USE_DEPENDENCY(cube-stacks "yaml-cpp")

TARGET_LINK_LIBRARIES(cube-stacks optimized manifolds debug manifolds_d)

TARGET_LINK_LIBRARIES(cube-stacks ${Boost_LIBRARIES})

SET(INSTALL_PATH include/cube-stacks)
INSTALL(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
INSTALL(TARGETS cube-stacks DESTINATION "lib")


