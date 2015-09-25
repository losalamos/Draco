# ----------------------------------------------------------------------
# Cray systems (Cielo, Cielito, Trinitite, Trinity)
#
# cmake -C ~/draco/config/CrayConfig.cmake [configure options] source_dir
# ----------------------------------------------------------------------

# Don't assume Linux to avoid -rdynamic flag when linking...
#set( CMAKE_SYSTEM_NAME Catamount CACHE STRING "description")

# Keyword for creating new libraries (STATIC or SHARED).
#set( DRACO_LIBRARY_TYPE STATIC CACHE STRING "description")

# Remove '-rdynamic' from link_flags
#foreach( lang CXX CC Fortran )
#  set(CMAKE_SHARED_LIBRARY_LINK_${lang}_FLAGS "" CACHE STRING "description" )
#endforeach()

# Use Cray provided compiler wrappers
set( CMAKE_C_COMPILER cc CACHE STRING "description" )
set( CMAKE_CXX_COMPILER CC CACHE STRING "description" )
set( CMAKE_Fortran_COMPILER ftn CACHE STRING "description" )

# Help find vendor software
set( MPIEXEC "aprun" CACHE STRING "description")
set( MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "description")

set( MPI_C_LIBRARIES "" CACHE STRING "description")
set( MPI_CXX_LIBRARIES "" CACHE STRING "description")
set( MPI_Fortran_LIBRARIES "" CACHE STRING "description")

set( MPI_C_INCLUDE_PATH "" CACHE STRING "description")
set( MPI_CXX_INCLUDE_PATH "" CACHE STRING "description")
set( MPI_Fortran_INCLUDE_PATH "" CACHE STRING "description")

set( MPI_CXX_COMPILER "" CACHE STRING "description")
set( MPI_C_COMPILER "" CACHE STRING "description")
set( MPI_Fortran_COMPILER "" CACHE STRING "description")

# Allinea MAP
# set( CMAKE_EXE_LINKER_FLAGS "" CACHE STRING "description" )