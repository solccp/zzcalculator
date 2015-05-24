# Project name is not mandatory, but you should use it
project(ZZCalculator Fortran)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)
 
# States that CMake required version must be >= 2.6
cmake_minimum_required(VERSION 2.6)

option(PROFILING "set profiling flags" OFF)
option(USE_OPENMP "enable OpenMP" OFF)

option(USE_MPI "Build parallel version with MPI library" OFF)

# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules (eg MacroOutOfSourceBuild.cmake):
set(CMAKE_MODULE_PATH ${ZZCalculator_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})
 
# The macro below forces the build directory to be different from source directory:
include(MacroOutOfSourceBuild)
 
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")

#find_package( Boost 1.34 COMPONENTS serialization REQUIRED )

add_subdirectory(src)
#SET_TARGET_PROPERTIES(ZZCalculator PROPERTIES LINKER_LANGUAGE Fortran)

enable_testing()
add_subdirectory(test)




SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_GENERATOR "STGZ;DEB")
SET(CPACK_PACKAGE_NAME "ZZCalculator")
SET(CPACK_PACKAGE_VENDOR "NCTU")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
#SET(CPACK_RESOURCE_FILE_README "README")
#SET(CPACK_SOURCE_GENERATOR "TGZ;TZ")
#SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_SYSTEM_NAME "Linux-amd64")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Chien-Pin Chou")
INCLUDE(CPack)
