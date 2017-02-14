INCLUDE_DIRECTORIES(BEFORE ${Boost_INCLUDE_DIR})

SET(SOURCES
  BoxTrajProblem.cc
  BoxesHullTrajProblem.cc
  utils/IndexManager.cc
  utils/ProblemConfig.cc
  utils/quat2mat.cc
  utils/Printer.cc
  functions/BoxAbovePlan.cc
  functions/BoxAboveFixedPlan.cc
  )
SET(HEADERS
  ../include/feet-trajectory/utils/Box.hh
  ../include/feet-trajectory/utils/FixedPlan.hh
  ../include/feet-trajectory/utils/Plan.hh
  ../include/feet-trajectory/utils/PlanForHull.hh
  ../include/feet-trajectory/utils/ProblemConfig.hh
  ../include/feet-trajectory/utils/ProblemInfo.hh
  ../include/feet-trajectory/utils/quat2mat.hh
  ../include/feet-trajectory/utils/Printer.hh
  ../include/feet-trajectory/functions/BoxAbovePlan.hh
  ../include/feet-trajectory/functions/BoxAboveFixedPlan.hh
  ../include/feet-trajectory/BoxTrajProblem.hh
  ../include/feet-trajectory/BoxesHullTrajProblem.hh
  )

ADD_LIBRARY(feet-trajectory SHARED ${SOURCES} ${HEADERS})

PKG_CONFIG_USE_DEPENDENCY(feet-trajectory roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(feet-trajectory roboptim-core-manifold)
PKG_CONFIG_USE_DEPENDENCY(feet-trajectory "yaml-cpp")

TARGET_LINK_LIBRARIES(feet-trajectory optimized manifolds debug manifolds_d)

TARGET_LINK_LIBRARIES(feet-trajectory ${Boost_LIBRARIES})

SET(INSTALL_PATH include/feet-trajectory)
INSTALL(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
INSTALL(TARGETS feet-trajectory DESTINATION "lib")


