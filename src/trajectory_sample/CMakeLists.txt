
add_library(trajectory_sample STATIC
  src/CartesianSample.cpp
  src/CurvilinearSample.cpp
  src/TrajectorySample.cpp
  src/FillCoordinates.cpp
  src/CoordinateSystemWrapper.cpp
)

# Find Eigen3
# find_package(Eigen3 REQUIRED NO_MODULE)

target_include_directories(trajectory_sample PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  # ${CMAKE_CURRENT_SOURCE_DIR}/../../include  
  # ${CMAKE_SOURCE_DIR}/eigen-3.4.0
  # ${CMAKE_SOURCE_DIR}/third_party/crdc/include
)

target_link_libraries(trajectory_sample PUBLIC crccosy cost_functions Eigen3::Eigen)
# target_include_directories(trajectory_sample PRIVATE SYSTEM ${EIGEN3_INCLUDE_DIR})
# target_include_directories(trajectory_sample INTERFACE 
#     ../../eigen-3.4.0
# )
