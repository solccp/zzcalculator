# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
 
# For a large number of source files you can create it in a simpler way
# using file() function:
file(GLOB ZZ_polynomial_SOURCES *.f90)

# Create an executable file called helloworld from sources:
add_executable(ZZ_polynomial ${ZZ_polynomial_SOURCES})

#set_target_properties(ZZ_polynomial PROPERTIES LINKER_LANGUAGE Fortran)
# test whether this is an intel compiler, then use -static-intel when linking
if( CMAKE_Fortran_COMPILER_ID STREQUAL Intel )
    if (CMAKE_BUILD_TYPE STREQUAL Debug)
        set(CMAKE_Fortran_FLAGS_DEBUG "")
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "-static-intel -static-libgcc")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O0 -g3 -trace -debug all -check all")
    endif (CMAKE_BUILD_TYPE STREQUAL Debug)
    
    if (CMAKE_BUILD_TYPE STREQUAL Release)
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "-static-intel -static-libgcc -ipo")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-ipo -O3 -funroll-loops")
        set(CMAKE_Fortran_FLAGS_RELEASE "")
    endif (CMAKE_BUILD_TYPE STREQUAL Release)
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
endif( CMAKE_Fortran_COMPILER_ID STREQUAL Intel )

if( CMAKE_Fortran_COMPILER_ID STREQUAL G95 )
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
    if (CMAKE_BUILD_TYPE STREQUAL Debug)
        set(CMAKE_Fortran_FLAGS_DEBUG "")
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O0 -g3")
    endif(CMAKE_BUILD_TYPE STREQUAL Debug)
    if (CMAKE_BUILD_TYPE STREQUAL Release)
        set(CMAKE_Fortran_FLAGS_RELEASE "")
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O3 -funroll-loops")
    endif(CMAKE_BUILD_TYPE STREQUAL Release)
endif( CMAKE_Fortran_COMPILER_ID STREQUAL G95 )
