# Project name is not mandatory, but you should use it
project(aromaticity Fortran)
 
# States that CMake required version must be >= 2.6
cmake_minimum_required(VERSION 2.6)

 
# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules (eg MacroOutOfSourceBuild.cmake):
set(CMAKE_MODULE_PATH ${aromaticity_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})
 
# The macro below forces the build directory to be different from source directory:
include(MacroOutOfSourceBuild)
 
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif (NOT CMAKE_BUILD_TYPE)


add_subdirectory(src)

enable_testing()
add_subdirectory(test)

#file(GLOB TEST_FILES test/*.in)

#foreach(loop_var IN LISTS TEST_FILES)
#    get_filename_component(FNAME ${loop_var} NAME_WE)
#    MESSAGE(${FNAME})
#    add_test(NAME ${FNAME}
#        COMMAND ${CMAKE_COMMAND}
#        -DTEST_PROG=$<TARGET_FILE:aromaticity>
#        -DTEST_NAME=${FNAME}
#        -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}/test
#        -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/script/runtest.cmake)

#endforeach()

