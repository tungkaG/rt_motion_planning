
add_library(initial_state STATIC
  src/InputData.cpp
)

target_include_directories(initial_state PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include  
)

target_include_directories(initial_state PRIVATE SYSTEM ${EIGEN3_INCLUDE_DIR})
