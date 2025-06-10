
add_library(trajectory_handler STATIC
  src/TrajectoryHandler.cpp
)

# find_package(Eigen3 REQUIRED NO_MODULE)

target_include_directories(trajectory_handler PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  # ${CMAKE_CURRENT_SOURCE_DIR}/../../include  
  # ${CMAKE_SOURCE_DIR}/eigen-3.4.0
)

target_link_libraries(trajectory_handler PUBLIC
  Eigen3::Eigen
  feasibility_functions
  trajectory_sample
  cost_functions
)
# target_include_directories(trajectory_sample PRIVATE SYSTEM ${EIGEN3_INCLUDE_DIR})
# target_include_directories(trajectory_sample INTERFACE 
#     ../../eigen-3.4.0
# )
