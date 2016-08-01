INCLUDE_DIRECTORIES(BEFORE ${Boost_INCLUDE_DIR})

SET(SOURCES
  CubeStackProblemOnManifold.cc
  CubeStackProblemOnSO3noS2.cc
  CubeStackProblemOnR.cc
  RoboptimCubeStackProblem.cc
  functions/Cost.cc
  functions/TotalCost.cc
  functions/CubeAboveFixedPlan.cc
  functions/CubeAboveFixedPlanRoboptim.cc
  functions/Norm1Quaternion.cc
  functions/Norm1Vector.cc
  functions/CubeAbovePlan.cc
  functions/CubeAbovePlanRoboptim.cc
  utils/ProblemConfig.cc
  utils/IndexManager.cc
  utils/quat2mat.cc
  )
SET(HEADERS
  ../include/cube-stacks/CubeStackProblemOnManifold.hh
  ../include/cube-stacks/CubeStackProblemOnSO3noS2.hh
  ../include/cube-stacks/CubeStackProblemOnR.hh
  ../include/cube-stacks/RoboptimCubeStackProblem.hh
  ../include/cube-stacks/functions/Cost.hh
  ../include/cube-stacks/functions/TotalCost.hh
  ../include/cube-stacks/functions/CubeAboveFixedPlan.hh
  ../include/cube-stacks/functions/CubeAboveFixedPlanRoboptim.hh
  ../include/cube-stacks/functions/Norm1Quaternion.hh
  ../include/cube-stacks/functions/Norm1Vector.hh
  ../include/cube-stacks/functions/CubeAbovePlan.hh
  ../include/cube-stacks/functions/CubeAbovePlanRoboptim.hh
  ../include/cube-stacks/utils/Cube.hh
  ../include/cube-stacks/utils/IndexManager.hh
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


