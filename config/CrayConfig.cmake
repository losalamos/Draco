# ----------------------------------------------------------------------
# Cray systems (Cielo, Cielito, Trinitite, Trinity)
#
# cmake -C ~/draco/config/CrayConfig.cmake [configure options] source_dir
# ----------------------------------------------------------------------

# Don't assume Linux to avoid -rdynamic flag when linking...
set( CMAKE_SYSTEM_NAME Catamount CACHE STRING "comp")

# Keyword for creating new libraries (STATIC or SHARED).
set( DRACO_LIBRARY_TYPE STATIC CACHE STRING "comp")

# Use Cray provided compiler wrappers
set( CMAKE_C_COMPILER cc CACHE STRING "comp" )
set( CMAKE_CXX_COMPILER CC CACHE STRING "comp" )
set( CMAKE_Fortran_COMPILER ftn CACHE STRING "comp" )

# Help find vendor software
set( MPIEXEC "aprun" CACHE STRING "comp")
set( MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "comp")

set( MPI_C_LIBRARIES "" CACHE STRING "comp")
set( MPI_CXX_LIBRARIES "" CACHE STRING "comp")
set( MPI_Fortran_LIBRARIES "" CACHE STRING "comp")

set( MPI_C_INCLUDE_PATH "" CACHE STRING "comp")
set( MPI_CXX_INCLUDE_PATH "" CACHE STRING "comp")
set( MPI_Fortran_INCLUDE_PATH "" CACHE STRING "comp")

set( MPI_CXX_COMPILER "" CACHE STRING "comp")
set( MPI_C_COMPILER "" CACHE STRING "comp")
set( MPI_Fortran_COMPILER "" CACHE STRING "comp")
