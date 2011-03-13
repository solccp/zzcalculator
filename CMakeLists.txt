# Project name is not mandatory, but you should use it
project(ZZ_polynomial Fortran)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)
 
# States that CMake required version must be >= 2.6
cmake_minimum_required(VERSION 2.6)

option(PROFILING "set profiling flags" OFF)
option(PARALLEL "set parallel flags" OFF)

# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules (eg MacroOutOfSourceBuild.cmake):
set(CMAKE_MODULE_PATH ${ZZ_polynomial_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})
 
# The macro below forces the build directory to be different from source directory:
include(MacroOutOfSourceBuild)
 
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")


add_subdirectory(src)

enable_testing()
add_subdirectory(test)


