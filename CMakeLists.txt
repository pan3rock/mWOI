  
cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0048 NEW)

project(mWOI VERSION 0.1)
set(CMAKE_CXX_STANDARD 17)

# set(CMAKE_CXX_FLAGS "-Wall -g -fPIC -O2 -march=native")
set(CMAKE_CXX_FLAGS "-Wall -g -fPIC")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# fmt
add_subdirectory(${CMAKE_SOURCE_DIR}/external/fmt EXCLUDE_FROM_ALL)

# yaml-cpp
add_subdirectory(${CMAKE_SOURCE_DIR}/external/yaml-cpp EXCLUDE_FROM_ALL)

# CLI11
add_subdirectory(${CMAKE_SOURCE_DIR}/external/CLI11 EXCLUDE_FROM_ALL)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(src)
add_subdirectory(test)

enable_testing()
add_test(NAME test_all COMMAND unit_test)