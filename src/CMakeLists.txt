
set(CMAKE_INCLUDE_CURRENT_DIR ON)
 
file(GLOB ZZ_polynomial_SOURCES *.f90 *.f95)

if(USE_MPI STREQUAL YES)
    LIST(APPEND ZZ_polynomial_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mpi/zhang.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mpi_global.f90)
    message("Build ZZ_polynomial with MPI support.")
else()
    #source_group("special" FILES non_mpi/*.f90)
    LIST(APPEND ZZ_polynomial_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/non_mpi/zhang.f90)
endif()


# Create an executable file called helloworld from sources:
find_package(LAPACK REQUIRED)

add_executable(ZZ_polynomial ${ZZ_polynomial_SOURCES})
target_link_libraries(ZZ_polynomial ${LAPACK_LIBRARIES})
if(USE_MPI STREQUAL YES)
    set_target_properties(ZZ_polynomial PROPERTIES OUTPUT_NAME "ZZ_polynomial_mpi")
elseif(USE_OPENMP STREQUAL YES)
    set_target_properties(ZZ_polynomial PROPERTIES OUTPUT_NAME "ZZ_polynomial_openmp")
endif()


#message ( "${CMAKE_Fortran_COMPILER_ID}"  )
#set_target_properties(ZZ_polynomial PROPERTIES LINKER_LANGUAGE Fortran)
# test whether this is an intel compiler, then use -static-intel when linking
#if( "${CMAKE_Fortran_COMPILER}" MATCHES "ifort" )
if( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
    if (CMAKE_BUILD_TYPE STREQUAL Debug)
        set(CMAKE_Fortran_FLAGS_DEBUG "")
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "-static-intel -static-libgcc ")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O0 -g3 -trace -debug all -check all")
    endif (CMAKE_BUILD_TYPE STREQUAL Debug)
    
    if (CMAKE_BUILD_TYPE STREQUAL Release)
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "-static-intel -static-libgcc -ipo ")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-ipo -O3 -funroll-loops -no-prec-div")
        if(USE_OPENMP STREQUAL YES)
            GET_TARGET_PROPERTY(COMFLAGS ZZ_polynomial COMPILE_FLAGS)
            GET_TARGET_PROPERTY(LINKFLAGS ZZ_polynomial LINK_FLAGS)
            set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "${LINKFLAGS} -openmp")
            set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "${COMFLAGS} -openmp")
        endif()
#        set(CMAKE_FLAGS_RELEASE "-ipo -O3 -funroll-loops -no-prec-div ")
        if (PROFILING STREQUAL YES)
            GET_TARGET_PROPERTY(COMFLAGS ZZ_polynomial COMPILE_FLAGS)
            GET_TARGET_PROPERTY(LINKFLAGS ZZ_polynomial LINK_FLAGS)
            set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "${COMFLAGS} -pg")
            set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "${COMFLAGS} -pg")
        endif()
        
    endif (CMAKE_BUILD_TYPE STREQUAL Release)
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
endif()

if (PROFILING)
    set (CMAKE_COMPILE_FLAGS "${CMAKE_COMPILE_FLAGS} -pg")
endif()
if( CMAKE_Fortran_COMPILER_ID STREQUAL G95 )
    set(CMAKE_SHARED_LIBRARY_LINK_FLAGS "")
    if (CMAKE_BUILD_TYPE STREQUAL Debug)
        set(CMAKE_FLAGS_DEBUG "")
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O0 -g3 -ftrace=full")
    endif(CMAKE_BUILD_TYPE STREQUAL Debug)
    if (CMAKE_BUILD_TYPE STREQUAL Release)
        set(CMAKE_FLAGS_RELEASE "")
        set_target_properties(ZZ_polynomial PROPERTIES LINK_FLAGS "")
        set_target_properties(ZZ_polynomial PROPERTIES COMPILE_FLAGS "-O3 -funroll-loops")
    endif(CMAKE_BUILD_TYPE STREQUAL Release)
endif( CMAKE_Fortran_COMPILER_ID STREQUAL G95 )




INSTALL(TARGETS ZZ_polynomial
  RUNTIME DESTINATION bin
)
