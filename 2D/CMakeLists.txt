cmake_minimum_required(VERSION 3.16)
project(Radial_migration VERSION 0.0.0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# You can set this flag to `-g -Wall` when developing to enable debugging
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

#Parallelization will be done later
# Find OpenMP and MPI
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(Catch2 3 REQUIRED)

include(CTest)

# Set include directory and add source directory
add_subdirectory(Log_pot)
add_subdirectory(VecUtils)
