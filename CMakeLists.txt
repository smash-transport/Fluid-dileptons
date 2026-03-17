cmake_minimum_required(VERSION 3.16)

project(DileptonHydroLib VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# sources
set(SOURCES
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

add_library(dilepton_hydro ${SOURCES} ${HEADERS})

target_include_directories(dilepton_hydro PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_features(dilepton_hydro PUBLIC cxx_std_17)

option(BUILD_EXAMPLE "Build example executable" ON)
if(BUILD_EXAMPLE AND EXISTS ${PROJECT_SOURCE_DIR}/example.cpp)
    add_executable(dilepton_example example.cpp)
    target_link_libraries(dilepton_example PRIVATE dilepton_hydro)
endif()