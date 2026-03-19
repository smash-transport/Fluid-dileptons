
cmake_minimum_required(VERSION 3.16)

project(TestDileptonHydro LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(FLUID_DILEPTONS_ROOT "" CACHE PATH "Path to fluid_dileptons library")

if(NOT FLUID_DILEPTONS_ROOT)
    message(FATAL_ERROR "Please set -DFLUID_DILEPTONS_ROOT=/path/to/library")
endif()

add_subdirectory(${FLUID_DILEPTONS_ROOT} ${CMAKE_BINARY_DIR}/fluid_dileptons_build)

add_executable(example_fluid_dileptons example.cpp)

target_link_libraries(example_fluid_dileptons PRIVATE fluid_dileptons)
