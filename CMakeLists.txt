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

# --- Dependencies ---

# OpenMP
find_package(OpenMP REQUIRED)

# HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSlib REQUIRED htslib)

# Eigen3
find_package(Eigen3 3.3 QUIET)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found via find_package, trying to find header path...")
    # Fallback for common locations
    if(EXISTS "/usr/include/eigen3")
        set(EIGEN3_INCLUDE_DIRS "/usr/include/eigen3")
        message(STATUS "Found Eigen3 at /usr/include/eigen3")
    else()
        message(WARNING "Eigen3 not found. Please install eigen3 or set EIGEN3_INCLUDE_DIRS.")
    endif()
else()
    set(EIGEN3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIRS})
endif()

# --- Includes ---
include_directories(include)
include_directories(${EIGEN3_INCLUDE_DIRS})

# --- Library ---
add_library(inter_sub_mod_core
    src/core/Config.cpp
    src/core/SomaticSnv.cpp
    src/utils/Logger.cpp
)

target_link_libraries(inter_sub_mod_core PUBLIC
    OpenMP::OpenMP_CXX
    ${HTSlib_LIBRARIES}
)
target_include_directories(inter_sub_mod_core PUBLIC 
    ${HTSlib_INCLUDE_DIRS}
)

# --- Main Executable ---
add_executable(inter_sub_mod src/main.cpp)
target_link_libraries(inter_sub_mod PRIVATE inter_sub_mod_core)

# --- Dev Tools / Ad-hoc Tests ---
add_executable(test_snv_loading src/test_snv_loading.cpp)
target_link_libraries(test_snv_loading PRIVATE inter_sub_mod_core)

# --- Tests ---
enable_testing()

# Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Check if we can fetch (assuming network might be an issue, we'll wrap in try/catch logic conceptually or just try)
# If offline, this will fail. I'll assume online access.
FetchContent_MakeAvailable(googletest)

add_executable(run_tests tests/test_config.cpp)
target_link_libraries(run_tests PRIVATE inter_sub_mod_core GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(run_tests)

