# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Project name and C++ standard
project(PMM_UAV_Planner LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set flags for compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

add_definitions(-DROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\")

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include /usr/include/eigen3)

# Source files
set(SOURCE_FILES
        src/pmm_mg_trajectory3d.cpp
        src/pmm_trajectory.cpp
        src/pmm_trajectory3d.cpp
        src/common.cpp
)

# Create executable
add_library(${PROJECT_NAME} ${SOURCE_FILES})

# Link libraries (pthread and yaml-cpp)
target_link_libraries(${PROJECT_NAME} pthread yaml-cpp)

add_executable(main
    src/main.cpp
)

target_link_libraries(main ${PROJECT_NAME} pthread yaml-cpp)

# Create an object directory for compiled files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/obj)

# Optional: If you have test files, you can add them here
# For example:
# set(TEST_FILES tests/test_pmm_mg_trajectory3d.cpp)
# add_executable(test_${PROJECT_NAME}
