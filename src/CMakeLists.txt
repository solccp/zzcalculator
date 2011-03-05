# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
 
# Create a variable called helloworld_SOURCES containing all .cpp files:
#set(aromaticity_SOURCES )
 
# For a large number of source files you can create it in a simpler way
# using file() function:
file(GLOB aromaticity_SOURCES *.f90)
 
# Create an executable file called helloworld from sources:
add_executable(aromaticity ${aromaticity_SOURCES})
# test whether this is an intel compiler, then use -static-intel when linking
if( CMAKE_Fortran_COMPILER_ID STREQUAL Intel )
#  set_target_properties(aromaticity PROPERTIES LINK_FLAGS "-static-intel -static-libgcc -ipo")
#  set_target_properties(aromaticity PROPERTIES COMPILE_FLAGS "-ipo -O3 -funroll-loops")
  #set(CMAKE_F_FLAGS_DEBUG -static-intel)
  set(CMAKE_F_FLAGS_DEBUG -static-libgcc)

  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
endif( CMAKE_Fortran_COMPILER_ID STREQUAL Intel )

if( CMAKE_Fortran_COMPILER_ID STREQUAL GNU )
  set_target_properties(aromaticity PROPERTIES LINK_FLAGS "-static-libgfortran -static-libgcc")
  set_target_properties(aromaticity PROPERTIES COMPILE_FLAGS "-O3 -g3 -funroll-loops")
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
endif( CMAKE_Fortran_COMPILER_ID STREQUAL GNU )
