# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
 
# For a large number of source files you can create it in a simpler way
# using file() function:
file(GLOB  ZZ_polynomial_SOURCES *.f90 *.cpp)
#source_group("common" FILES *.f90)

if(USE_MPI STREQUAL ON)
    #source_group("special" FILES mpi/*.f90)
    LIST(APPEND ZZ_polynomial_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mpi/zhang.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mpi_global.f90)
    message("Build ZZ_polynomial with MPI support.")
else()
    #source_group("special" FILES non_mpi/*.f90)
    LIST(APPEND ZZ_polynomial_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/non_mpi/zhang.f90)
endif()


# Create an executable file called helloworld from sources:
add_executable(ZZ_polynomial ${ZZ_polynomial_SOURCES})
if(USE_MPI STREQUAL ON)
    set_target_properties(ZZ_polynomial PROPERTIES OUTPUT_NAME "ZZ_polynomial_mpi")
endif()

#if (PROFILING)
#    set (CMAKE_COMPILE_FLAGS "-pg")
#endif()

#set_target_properties(ZZ_polynomial PROPERTIES LINKER_LANGUAGE Fortran)
# test whether this is an intel compiler, then use -static-intel when linking
if( CMAKE_Fortran_COMPILER_ID STREQUAL Intel )
    if (CMAKE_BUILD_TYPE STREQUAL Debug)
        set(CMAKE_Fortran_FLAGS_DEBUG "")
        set_target_properties(ZZ_polynomial PROPERTIES Fortran_LINK_FLAGS "-static-intel -static-libgcc")
        set_target_properties(ZZ_polynomial PROPERTIES Fortran_COMPILE_FLAGS "-O0 -g3 -trace -debug all -check all")
    endif (CMAKE_BUILD_TYPE STREQUAL Debug)
    
    if (CMAKE_BUILD_TYPE STREQUAL Release)
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "-static-intel -static-libgcc -ipo")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O3 ")
#        set_target_properties(ZZ_polynomial PROPERTIES Fortran_COMPILE_FLAGS "-ipo -O3 -funroll-loops -no-prec-div")
        set(CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -funroll-loops -no-prec-div ")
        if (PROFILING)
            set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O3 -pg")
            set(CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -funroll-loops -no-prec-div -pg")
        endif()
        
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
        set_target_properties(ZZ_polynomial PROPERTIES Fortran_COMPILE_FLAGS "-O3 -funroll-loops")
    endif(CMAKE_BUILD_TYPE STREQUAL Release)
endif( CMAKE_Fortran_COMPILER_ID STREQUAL G95 )


target_link_libraries(ZZ_polynomial
  ${Boost_SERIALIZATION_LIBRARY}
)

INSTALL(TARGETS ZZ_polynomial
  RUNTIME DESTINATION bin
)
