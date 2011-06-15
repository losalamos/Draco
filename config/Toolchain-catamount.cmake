# -------- Toolchain-catamount.cmake ---------

# Ref: http://www.vtk.org/Wiki/CMake_Cross_Compiling
#      http://www.vtk.org/Wiki/CmakeCrayXt3

# Specify the target system (this allows cross-compiling)
SET(CMAKE_SYSTEM_NAME Catamount)

# specify the cross compiler
SET(CMAKE_C_COMPILER   cc) # --target=catamount)
SET(CMAKE_CXX_COMPILER CC) # --target=catamount)

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /opt/xt-pe/default
    /opt/xt-mpt/default/mpich2-64/GP2
    $ENV{MPICH_DIR}/lib
# DRACO_DIR?
# VENDOR_DIR?    
  )

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment.  Can be set to NEVER, ONLY or BOTH (default)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)

# Help find vendor software
set( MPI_INC_DIR $ENV{MPICH_DIR}/include )
set( MPI_LIB_DIR $ENV{MPICH_DIR}/lib )
set( DRACO_LIBRARY_TYPE "STATIC" )
set( MPIEXEC "/usr/bin/aprun" )
set( MPIEXEC_NUMPROC_FLAG "-n")