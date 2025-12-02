cmake_minimum_required(VERSION 3.14)
project(InterSubMod CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Compiler Options ---
add_compile_options(-Wall -Wextra -pedantic)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0)
    add_definitions(-DDEBUG)
else()
    add_compile_options(-O3 -march=native)
    add_definitions(-DNDEBUG)
endif()

# --- Output Directory ---
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Jemalloc Integration ---
include(ExternalProject)
set(JEMALLOC_PREFIX ${CMAKE_BINARY_DIR}/external/jemalloc)
set(JEMALLOC_LIB_DIR ${JEMALLOC_PREFIX}/lib)
set(JEMALLOC_INCLUDE_DIR ${JEMALLOC_PREFIX}/include)

# Create directories to satisfy CMake checks
file(MAKE_DIRECTORY ${JEMALLOC_INCLUDE_DIR})
file(MAKE_DIRECTORY ${JEMALLOC_LIB_DIR})

# Download and build jemalloc statically
ExternalProject_Add(jemalloc_external
    URL https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2
    PREFIX ${JEMALLOC_PREFIX}
    # Configure: disable shared, enable static, no prefix (replace system malloc)
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-shared --enable-static --with-jemalloc-prefix= 
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)

# Define imported library for jemalloc
add_library(jemalloc STATIC IMPORTED)
set_target_properties(jemalloc PROPERTIES
    IMPORTED_LOCATION ${JEMALLOC_LIB_DIR}/libjemalloc.a
    INTERFACE_INCLUDE_DIRECTORIES ${JEMALLOC_INCLUDE_DIR}
)
add_dependencies(jemalloc jemalloc_external)

# Enable jemalloc usage in code
add_definitions(-DUSE_JEMALLOC)

# --- Dependencies ---

# OpenMP
find_package(OpenMP REQUIRED)

# HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSlib REQUIRED htslib)

# Eigen3
include(FetchContent)
FetchContent_Declare(
  eigen3
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
)
FetchContent_MakeAvailable(eigen3)
include_directories(${eigen3_SOURCE_DIR})

# --- Includes ---
include_directories(include)
include_directories(${EIGEN3_INCLUDE_DIRS})
# jemalloc include dir will be added via target usage, but global include helps too
include_directories(${JEMALLOC_INCLUDE_DIR}) 

# --- Library ---
add_library(inter_sub_mod_core
    src/core/Config.cpp
    src/core/SomaticSnv.cpp
    src/core/BamReader.cpp
    src/core/ReadParser.cpp
    src/core/MethylationParser.cpp
    src/core/MatrixBuilder.cpp
    src/core/RegionProcessor.cpp
    src/core/DistanceMatrix.cpp
    src/utils/Logger.cpp
    src/utils/FastaReader.cpp
    src/io/RegionWriter.cpp
)

target_link_libraries(inter_sub_mod_core PUBLIC
    OpenMP::OpenMP_CXX
    ${HTSlib_LIBRARIES}
    jemalloc
    dl 
    pthread
)
target_include_directories(inter_sub_mod_core PUBLIC 
    ${HTSlib_INCLUDE_DIRS}
    ${JEMALLOC_INCLUDE_DIR}
)

# --- Main Executable ---
add_executable(inter_sub_mod src/main.cpp)
target_link_libraries(inter_sub_mod PRIVATE inter_sub_mod_core)

# --- Dev Tools / Ad-hoc Tests ---
add_executable(test_phase1_2 src/test/test_phase1_2.cpp)
target_link_libraries(test_phase1_2 PRIVATE inter_sub_mod_core)

add_executable(test_phase3 src/test/test_phase3.cpp)
target_link_libraries(test_phase3 PRIVATE inter_sub_mod_core)

add_executable(test_phase4_5 src/test/test_phase4_5.cpp)
target_link_libraries(test_phase4_5 PRIVATE inter_sub_mod_core)

# --- Tests ---
enable_testing()

# Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Run tests with custom main to track resources
add_executable(run_tests 
    tests/main_test.cpp 
    tests/test_config.cpp
    tests/test_bam_reader.cpp
    tests/test_distance_matrix.cpp
)
target_link_libraries(run_tests PRIVATE inter_sub_mod_core GTest::gtest)

include(GoogleTest)
gtest_discover_tests(run_tests)
