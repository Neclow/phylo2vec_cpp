cmake_minimum_required(VERSION 3.22.1)
project(phylo2vec)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Google Test using FetchContent
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.11.0 # dpkg -l libgtest-dev
)
FetchContent_MakeAvailable(googletest)

# Include cxxopts using FetchContent
FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG v3.1.1
)
FetchContent_MakeAvailable(cxxopts)

# # Find OpenMP package
# find_package(OpenMP REQUIRED)

# # Set compiler flags for OpenMP
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Source files
set(SOURCES
    src/phylo2vec.cpp
    src/main.cpp
)

# Test
set(TEST_SOURCES
    src/phylo2vec.cpp
    test/phylo2vec_test.cpp
)

# Main executable
add_executable(phylo2vec ${SOURCES})

# Include directories
target_include_directories(phylo2vec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Link against cxxopts
target_link_libraries(phylo2vec PRIVATE cxxopts::cxxopts) ## OpenMP::OpenMP_CXX)

# Test executable
add_executable(phylo2vec_test ${TEST_SOURCES})

# Link against Google Test and Google Mock
target_link_libraries(phylo2vec_test PRIVATE gtest_main)

# Optionally, add a test target (for running tests using CTest)
# include(CTest)

# Add tests
add_test(NAME phylo2vec_test COMMAND phylo2vec_test)
