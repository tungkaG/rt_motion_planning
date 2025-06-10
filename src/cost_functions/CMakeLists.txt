add_library(cost_functions STATIC
  src/CalculateAccelerationCost.cpp
  src/CalculateDistanceToReferencePathCost.cpp
  src/CalculateJerkCost.cpp
  src/CalculateLateralAccelerationCost.cpp
  src/CalculateLateralVelocityCost.cpp
  src/CalculateLongitudinalAccelerationCost.cpp
  src/CalculateLongitudinalVelocityCost.cpp
  src/CalculateNegativeAccelerationCost.cpp
  src/CalculateNegativeOrientationOffsetCost.cpp
  src/CalculateNegativeVelocityOffsetCost.cpp
  src/CalculateOrientationOffsetCost.cpp
  src/CalculatePositiveAccelerationCost.cpp
  src/CalculatePositiveOrientationOffsetCost.cpp
  src/CalculatePositiveVelocityOffsetCost.cpp
  src/CalculateVelocityOffsetCost.cpp
  src/util.cpp
)

# find_package(Eigen3 REQUIRED NO_MODULE)

target_include_directories(cost_functions PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  # ${CMAKE_SOURCE_DIR}/eigen-3.4.0
)

# target_include_directories(cost_functions PRIVATE SYSTEM ${EIGEN3_INCLUDE_DIR})

target_link_libraries(cost_functions PUBLIC Eigen3::Eigen trajectory_sample)
