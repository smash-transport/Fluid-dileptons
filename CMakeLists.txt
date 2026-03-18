cmake_minimum_required(VERSION 3.16)

project(FluidDileptonsLib VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# sources
set(SOURCES
    src/acceptance.cpp
    src/dilepton.cpp
    src/fourpion_mix.cpp
    src/output.cpp
    src/hydrocell.cpp
    src/setup.cpp
)

# Headers
set(HEADERS
    include/acceptance.h
    include/dilepton.h
    include/hydrocell.h
    include/output.h
    include/setup.h
    include/vector.h
)

add_library(fluid_dileptons ${SOURCES} ${HEADERS})

target_include_directories(fluid_dileptons PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_features(fluid_dileptons PUBLIC cxx_std_17)

option(BUILD_EXAMPLE "Build example executable" ON)
if(BUILD_EXAMPLE AND EXISTS ${PROJECT_SOURCE_DIR}/tools/example.cpp)
    add_executable(dilepton_example tools/example.cpp)
    target_link_libraries(dilepton_example PRIVATE fluid_dileptons)
endif()