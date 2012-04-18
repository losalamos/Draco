#-----------------------------*-cmake-*----------------------------------------#
# file   config/vendor_libraries.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 6
# brief  Setup Vendors
# note   Copyright (C) 2010-2012 LANS, LLC  
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#

#
# Look for any libraries which are required at the toplevel.
# 

include( FeatureSummary )

#------------------------------------------------------------------------------#
# Setup MPI when on Linux
#------------------------------------------------------------------------------#
macro( setupMPILibrariesUnix )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )
      # Try to find MPI in the default locations (look for mpic++ in PATH)
      # This module will set the following variables:
      #   MPI_FOUND                  TRUE have we found MPI
      #   MPI_COMPILE_FLAGS          Compilation flags for MPI programs
      #   MPI_INCLUDE_PATH           Include path(s) for MPI header
      #   MPI_LINK_FLAGS             Linking flags for MPI programs
      #   MPI_LIBRARY                First MPI library to link against (cached)
      #   MPI_EXTRA_LIBRARY          Extra MPI libraries to link against (cached)
      #   MPI_LIBRARIES              All libraries to link MPI programs against
      #   MPIEXEC                    Executable for running MPI programs
      #   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving it the
      #                              number of processors to run on
      #   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly before the
      #                              executable to run.
      #   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after all other flags.  
      
      # message(STATUS "Looking for MPI...")

      # Preserve data that may already be set.
      if( DEFINED MPIEXEC )
         set( MPIEXEC ${MPIEXEC} CACHE FILEPATH "Executable for running MPI programs." )
      endif()
      if( DEFINED MPIEXEC_MAX_NUMPROCS )
         set( MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS} CACHE STRING
            "Maximum number of processors available to run MPI applications.")
      endif()
      if( DEFINED MPIEXEC_NUMPROC_FLAG )
         set( MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG} CACHE
      STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes." )
      endif()

      # First attempt to find mpi
      find_package( MPI QUIET )

      # Second chance using $MPIRUN (old Draco setup format -- ask JDD).
      if( NOT ${MPI_FOUND} AND EXISTS "${MPIRUN}" )
#         message( STATUS "2nd attempt to find MPI: examine $MPIRUN" )
         set( MPIEXEC $ENV{MPIRUN} )
         find_package( MPI QUIET )
      endif()

      # Third chance using $MPI_INC_DIR and $MPI_LIB_DIR
      if( NOT ${MPI_FOUND} AND EXISTS "${MPI_LIB_DIR}" AND 
            EXISTS "${MPI_INC_DIR}" )
#         message( STATUS "3rd attempt to find MPI: examine $MPI_INC_DIR" )
         if( EXISTS "$ENV{MPI_INC_DIR}" AND "${MPI_INC_DIR}x" MATCHES "x" )
            set( MPI_INC_DIR $ENV{MPI_INC_DIR} )
         endif()
         if( EXISTS "$ENV{MPI_LIB_DIR}" AND "${MPI_LIB_DIR}x" MATCHES "x" )
            set( MPI_LIB_DIR $ENV{MPI_LIB_DIR} )
         endif()
         set( MPI_INCLUDE_PATH ${MPI_INC_DIR} )
         find_library( MPI_LIBRARY
            NAMES mpi mpich msmpi
            PATHS ${MPI_LIB_DIR} 
            ${MPICH_DIR}/lib
            )
         set( extra_libs mpi++ libopen-rte libopen-pal)
         unset( MPI_EXTRA_LIBRARY )
         foreach( lib ${extra_libs} )
            find_library( mpi_extra_lib_${lib}
               NAMES ${lib}
               HINTS ${MPI_LIB_DIR} 
               ${MPICH_DIR}/lib )
            mark_as_advanced( mpi_extra_lib_${lib} )
            if( EXISTS "${mpi_extra_lib_${lib}}" )
               list( APPEND MPI_EXTRA_LIBRARY ${tmp} )
            endif()
         endforeach()
         find_package( MPI QUIET )
         if( ${MPI_EXTRA_LIBRARY} MATCHES "NOTFOUND" )
            # do nothing
         else()
            list( APPEND MPI_LIBRAIES ${MPI_EXTRA_LIBRARY} )
         endif()
      endif()

      # Set Draco build system variables based on what we know about MPI.
      if( MPI_FOUND )
#         message( STATUS "MPI FOUND: Setting DRACO_C4 = MPI" )
         set( DRACO_C4 "MPI" )  
         set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMPI_SKIP_MPICXX" )
         if( NOT MPIEXEC )
            if( "${SITE}" MATCHES "c[it]" )
               set( MPIEXEC aprun )
            else()
               message( FATAL_ERROR 
                  "MPI found but mpirun not in PATH. Aborting" )
            endif()
         endif()
         # Try to find the fortran mpi library
         if( EXISTS ${MPI_LIB_DIR} )
            find_library( MPI_Fortran_LIB mpi_f77 HINTS ${MPI_LIB_DIR} )
            mark_as_advanced( MPI_Fortran_LIB )
         endif()
         
      else()
#         message( STATUS "MPI NOT FOUND: Setting DRACO_C4 = SCALAR" )
         set( DRACO_C4 "SCALAR" )
         set( MPI_INCLUDE_PATH "" CACHE FILEPATH "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARY "" CACHE FILEPATH "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARIES "" CACHE FILEPATH "no mpi library for scalar build." FORCE )
      endif()

      # Save the result in the cache file.
      set( DRACO_C4 "${DRACO_C4}" CACHE STRING 
         "C4 communication mode (SCALAR or MPI)" )
      # Provide a constrained pull down list in cmake-gui
      set_property( CACHE DRACO_C4 PROPERTY STRINGS SCALAR MPI )
      if( "${DRACO_C4}" STREQUAL "MPI"    OR 
            "${DRACO_C4}" STREQUAL "SCALAR" )
      else()
         message( FATAL_ERROR "DRACO_C4 must be either MPI or SCALAR" )
      endif()

      # Check flavor and add optional flags
      if( "${MPIEXEC}" MATCHES openmpi )
         set( MPI_FLAVOR "openmpi" CACHE STRING "Flavor of MPI." )
         execute_process(
            COMMAND ${MPIEXEC} --version
            ERROR_VARIABLE  DBS_MPI_VER
            )
         string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\1"
            DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
         string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\2"
            DBS_MPI_VER_MINOR ${DBS_MPI_VER} )

         # Ref: http://www.open-mpi.org/faq/?category=tuning#using-paffinity-v1.2
         # This is required on Turning when running 'ctest -j16'.  See
         # notes in component_macros.cmake.
        
         # --bind-to-core added in OpenMPI-1.4
         if( ${DBS_MPI_VER_MINOR} GREATER 3 )
            set( MPIEXEC_POSTFLAGS --mca mpi_paffinity_alone 0 CACHE
               STRING "extra mpirun flags (list)." FORCE)
            set( MPIEXEC_POSTFLAGS_STRING "--mca mpi_paffinity_alone 0" CACHE
               STRING "extra mpirun flags (string)." FORCE)
            # set( MPIEXEC_POSTFLAGS --bind-to-none --bycore CACHE
            #    STRING "extra mpirun flags (list)." FORCE)
            # set( MPIEXEC_POSTFLAGS_STRING "--bind-to-none --bycore" CACHE
            #    STRING "extra mpirun flags (string)." FORCE)
         else()
            set( MPIEXEC_POSTFLAGS --mca mpi_paffinity_alone 0 CACHE
               STRING "extra mpirun flags (list)." FORCE)
            set( MPIEXEC_POSTFLAGS_STRING "--mca mpi_paffinity_alone 0" CACHE
               STRING "extra mpirun flags (string)." FORCE)
         endif()
         mark_as_advanced( MPI_FLAVOR MPIEXEC_POSTFLAGS_STRING )
      endif()

      # Mark some of the variables created by the above logic as
      # 'advanced' so that they do not show up in the 'simple' ccmake 
      # view. 
      mark_as_advanced( MPI_EXTRA_LIBRARY MPI_LIBRARY )
      set( file_cmd ${file_cmd} CACHE INTERNAL "file command" )

   endif( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

endmacro()

#------------------------------------------------------------------------------
# Helper macros for LAPACK/Unix
#
# This module sets the following variables:
# LAPACK_FOUND - set to true if a library implementing the LAPACK
#         interface is found 
# LAPACK_LINKER_FLAGS - uncached list of required linker flags
#         (excluding -l and -L). 
# LAPACK_LIBRARIES - uncached list of libraries (using full path name)
#         to link against to use LAPACK 
#------------------------------------------------------------------------------
macro( setupLAPACKLibrariesUnix )
   if( NOT EXISTS ${LAPACK_lapack_LIBRARY} )

      set( BLA_STATIC ON )
      set( BLAS_REQUIRED "" )
      # set( ENV{BLA_VENDOR} "Generic" )
      
      if( EXISTS $ENV{LAPACK_LIB_DIR} )
         set( ENV{LD_LIBRARY_PATH}
            "$ENV{LAPACK_LIB_DIR}:$ENV{LD_LIBRARY_PATH}")
      endif()
      
      if( ${SITE} MATCHES "frost" )
         # This machine uses Intel MKL instead of LAPACK. 
         # See http://software.intel.com/sites/products/documentation/
         # hpc/compilerpro/en-us/fortran/lin/mkl/userguide.pdf
         set( ENV{BLA_VENDOR} "Intel10_64lp_gf_sequential" )
         find_package( BLAS ${BLAS_REQUIRED} QUIET )
         # We don't need mkl_lapack so just set LAPACK_LIBRARIES to
         # BLAS_LIBRARIES 
         set( LAPACK_LIBRARIES ${BLAS_LIBRARIES} )
         set( LAPACK_FOUND ON )
      else()
         if( EXISTS $ENV{LAPACK_LIB_DIR}/libblas.a )
            # avoid picking /usr/lib64/libblas.a
            set( BLAS_blas_LIBRARY $ENV{LAPACK_LIB_DIR}/libblas.a )
         endif()
         if( EXISTS $ENV{LAPACK_LIB_DIR}/liblapack.a )
            # avoid picking /usr/lib64/liblapack.a
            set( LAPACK_lapack_LIBRARY $ENV{LAPACK_LIB_DIR}/liblapack.a )
         endif()
         find_package( LAPACK ) 
      endif()
      
      if( LAPACK_FOUND )
         set( LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING 
            "lapack libs" )
         mark_as_advanced( LAPACK_LIBRARIES )
      endif()
      
   endif( NOT EXISTS ${LAPACK_lapack_LIBRARY} )

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesUnix )

   # GSL ----------------------------------------------------------------------
   # message( STATUS "Looking for GSL...")
   find_package( GSL QUIET )

   # GRACE ------------------------------------------------------------------
   find_package( Grace QUIET )
   set_package_properties( Grace PROPERTIES
      DESCRIPTION "A WYSIWYG 2D plotting tool."
      TYPE OPTIONAL
      PURPOSE "Required for bulding the plot2D component."
      )

   # CUDA ------------------------------------------------------------------
   find_package( CUDA QUIET )
   set_package_properties( CUDA PROPERTIES
      DESCRIPTION "Toolkit providing tools and libraries needed for GPU applications."
      TYPE OPTIONAL
      PURPOSE "Required for bulding a GPU enabled application." )
   if( NOT EXISTS ${CUDA_NVCC_EXECUTABLE} )
      set( CUDA_FOUND 0 )
   endif()
   if( CUDA_FOUND )
      set( HAVE_CUDA 1 )
      option( USE_CUDA "If CUDA is available, should we use it?" ON )
      set( CUDA_PROPAGATE_HOST_FLAGS OFF CACHE BOOL "blah" FORCE)
      set( CUDA_NVCC_FLAGS "-arch=sm_21" )
      set( cudalibs ${CUDA_CUDART_LIBRARY} )
   endif()
   mark_as_advanced( 
      CUDA_SDK_ROOT_DIR 
      CUDA_VERBOSE_BUILD
      CUDA_TOOLKIT_ROOT_DIR 
      CUDA_BUILD_CUBIN
      CUDA_BUILD_EMULATION
      )
endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesWindows )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )
      # Try to find MPI in the default locations (look for mpic++ in PATH)
      # This module will set the following variables:
      #   MPI_FOUND                  TRUE have we found MPI
      #   MPI_COMPILE_FLAGS          Compilation flags for MPI programs
      #   MPI_INCLUDE_PATH           Include path(s) for MPI header
      #   MPI_LINK_FLAGS             Linking flags for MPI programs
      #   MPI_LIBRARY                First MPI library to link against (cached)
      #   MPI_EXTRA_LIBRARY          Extra MPI libraries to link against (cached)
      #   MPI_LIBRARIES              All libraries to link MPI programs against
      #   MPIEXEC                    Executable for running MPI programs
      #   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving it the
      #                              number of processors to run on
      #   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly before the
      #                              executable to run.
      #   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after all other flags.  
      
      message(STATUS "Looking for MPI...")
      find_package( MPI )

      # Second chance using $MPIRUN (old Draco setup format -- ask JDD).
      if( NOT ${MPI_FOUND} AND EXISTS "${MPIRUN}" )
         set( MPIEXEC $ENV{MPIRUN} )
         find_package( MPI )
      endif()

      # Third chance using $MPI_INC_DIR and $MPI_LIB_DIR
      if( NOT ${MPI_FOUND} AND EXISTS "${MPI_LIB_DIR}" AND 
            EXISTS "${MPI_INC_DIR}" )
         if( EXISTS "$ENV{MPI_INC_DIR}" AND "${MPI_INC_DIR}x" MATCHES "x" )
            set( MPI_INC_DIR $ENV{MPI_INC_DIR} )
         endif()
         if( EXISTS "$ENV{MPI_LIB_DIR}" AND "${MPI_LIB_DIR}x" MATCHES "x" )
            set( MPI_LIB_DIR $ENV{MPI_LIB_DIR} )
         endif()
         set( MPI_INCLUDE_PATH ${MPI_INC_DIR} )
         find_library( MPI_LIBRARY
            NAMES mpi mpich msmpi
            PATHS ${MPI_LIB_DIR} 
            ${MPICH_DIR}/lib
            )
         set( extra_libs mpi++ libopen-rte libopen-pal)
         set( MPI_EXTRA_LIBRARY )
         foreach( lib ${extra_libs} )
            find_library( mpi_extra_lib_${lib}
               NAMES ${lib}
               HINTS ${MPI_LIB_DIR} 
               ${MPICH_DIR}/lib )
            mark_as_advanced( mpi_extra_lib_${lib} )
            if( EXISTS "${mpi_extra_lib_${lib}}" )
               list( APPEND MPI_EXTRA_LIBRARY ${tmp} )
            endif()
         endforeach()
         find_package( MPI )
         if( ${MPI_EXTRA_LIBRARY} MATCHES "NOTFOUND" )
            # do nothing
         else()
            list( APPEND MPI_LIBRAIES ${MPI_EXTRA_LIBRARY} )
         endif()
      endif()

      # Set Draco build system variables based on what we know about MPI.
      if( MPI_FOUND )
         set( DRACO_C4 "MPI" )  
         set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMPI_SKIP_MPICXX" )
         if( NOT MPIEXEC )
            message( FATAL_ERROR 
               "MPI found but mpirun not in PATH. Aborting" )
         endif()
         # Try to find the fortran mpi library
         if( EXISTS ${MPI_LIB_DIR} )
            find_library( MPI_Fortran_LIB mpi_f77 HINTS ${MPI_LIB_DIR} )
            mark_as_advanced( MPI_Fortran_LIB )
         endif()
         unset( C4_SCALAR )
      else()
         set( DRACO_C4 "SCALAR" )
         unset( C4_MPI )
         set( MPI_INCLUDE_PATH "" CACHE FILEPATH "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARY "" CACHE FILEPATH "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARIES "" CACHE FILEPATH "no mpi library for scalar build." FORCE )
      endif()

      # Save the result in the cache file.
      set( DRACO_C4 "${DRACO_C4}" CACHE STRING 
         "C4 communication mode (SCALAR or MPI)" )
      if( "${DRACO_C4}" STREQUAL "MPI"    OR 
            "${DRACO_C4}" STREQUAL "SCALAR" )
      else()
         message( FATAL_ERROR "DRACO_C4 must be either MPI or SCALAR" )
      endif()

      # Check flavor and add optional flags
      if( "${MPIEXEC}" MATCHES openmpi )
         set( MPI_FLAVOR "openmpi" CACHE STRING "Flavor of MPI." )
         # Ref: http://www.open-mpi.org/faq/?category=tuning#using-paffinity-v1.2
         # This is required on Turning when running 'ctest -j16'.  See
         # notes in component_macros.cmake.
         set( MPIEXEC_POSTFLAGS --mca mpi_paffinity_alone 0 CACHE
            STRING "extra mpirun flags (list)." FORCE)
         set( MPIEXEC_POSTFLAGS_STRING "--mca mpi_paffinity_alone 0" CACHE
            STRING "extra mpirun flags (string)." FORCE)
         mark_as_advanced( MPI_FLAVOR MPIEXEC_POSTFLAGS_STRING )
      endif()
      
   endif( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

   # LAPACK ------------------------------------------------------------------
   if( NOT EXISTS ${LAPACK_lapack_LIBRARY} )
      message( STATUS "Looking for LAPACK...")
      
      # Set the BLAS/LAPACK VENDOR.  
      # set( BLA_VENDOR "Generic" )
      
      # This module sets the following variables:
      # LAPACK_FOUND - set to true if a library implementing the LAPACK
      #              interface is found
      # LAPACK_LINKER_FLAGS - uncached list of required linker flags
      #              (excluding -l and -L).
      # LAPACK_LIBRARIES - uncached list of libraries (using full path
      #              name) to link against to use LAPACK

      # Use BLA_VENDOR and BLA_STATIC values from the BLAS section above.

      find_package( LAPACK ) # QUIET

      if( LAPACK_FOUND )
         set( LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING "lapack libs" )
         mark_as_advanced( LAPACK_LIBRARIES )
         message( STATUS "Found LAPACK: ${LAPACK_LIBRARIES}" )
      endif()
      
   endif( NOT EXISTS ${LAPACK_lapack_LIBRARY} )
   
   # GSL ---------------------------------------------------------------------
   message( STATUS "Looking for GSL...")
   set( GSL_INC_DIR "${VENDOR_DIR}/gsl/include" )
   set( GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib" )
   
   # Use static BLAS libraries
   set(GSL_STATIC ON)
   
   find_package( GSL REQUIRED )

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
# Assign here the library version to be used.
#------------------------------------------------------------------------------
macro( setVendorVersionDefaults )
    #Set the preferred search directories(ROOT)

    #Check that VENDOR_DIR is defined as a cache variable or as an
    #environment variable. If defined as both then take the
    #environment variable.

    # See if VENDOR_DIR is set.  Try some defaults if it is not set.
    if( NOT EXISTS "${VENDOR_DIR}" AND IS_DIRECTORY "$ENV{VENDOR_DIR}" )
        set( VENDOR_DIR $ENV{VENDOR_DIR} )
    endif()
    # If needed, try some obvious places.
    if( NOT EXISTS "${VENDOR_DIR}" )
        if( IS_DIRECTORY "/ccs/codes/radtran/vendors/Linux64" )
            set( VENDOR_DIR "/ccs/codes/radtran/vendors/Linux64" )
        endif()
        if( IS_DIRECTORY /usr/projects/draco/vendors )
            set( VENDOR_DIR /usr/projects/draco/vendors )
        endif()
        if( IS_DIRECTORY c:/vendors )
            set( VENDOR_DIR c:/vendors )
        elseif( IS_DIRECTORY c:/a/vendors )
            set( VENDOR_DIR c:/a/vendors )
        endif()
    endif()
    # Cache the result
    if( IS_DIRECTORY "${VENDOR_DIR}")
        set( VENDOR_DIR $ENV{VENDOR_DIR} CACHE PATH
          "Root directory where CCS-2 3rd party libraries are located." )
    else( IS_DIRECTORY "${VENDOR_DIR}")
        message( "
WARNING: VENDOR_DIR not defined locally or in user environment,
individual vendor directories should be defined." )
    endif()

    # Import environment variables related to vendors
    # 1. Use command line variables (-DLAPACK_LIB_DIR=<path>
    # 2. Use environment variables ($ENV{LAPACK_LIB_DIR}=<path>)
    # 3. Try to find vendor in $VENDOR_DIR
    # 4. Don't set anything and let the user set a value in the cache
    #    after failed 1st configure attempt.
    if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY $ENV{LAPACK_LIB_DIR} )
        set( LAPACK_LIB_DIR $ENV{LAPACK_LIB_DIR} )
        set( LAPACK_INC_DIR $ENV{LAPACK_INC_DIR} )
    endif()
    if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/LAPACK/lib )
        set( LAPACK_LIB_DIR "${VENDOR_DIR}/LAPACK/lib" )
        set( LAPACK_INC_DIR "${VENDOR_DIR}/LAPACK/include" )
        if( WIN32 )
            # cleanup LIB
            unset(newlibpath)
            foreach( path $ENV{LIB} )
                if( newlibpath )
                    set( newlibpath "${newlibpath};${path}" )
                else()
                    set( newlibpath "${path}" )
                endif()
            endforeach()
            set(haslapackpath FALSE)
            foreach( path ${newlibpath} )
                if( "${LAPACK_LIB_DIR}" STREQUAL "${path}" )
                    set( haslapackpath TRUE ) 
                endif()
            endforeach()
            if( NOT ${haslapackpath} )
                set( newlibpath "${newlibpath};${LAPACK_LIB_DIR}" )
                set( ENV{LIB} "${newlibpath}" )
            endif()
            #message("LIB = $ENV{LIB}")
        endif()
    endif()
    if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/clapack/lib )
        set( LAPACK_LIB_DIR "${VENDOR_DIR}/clapack/lib" )
        set( LAPACK_INC_DIR "${VENDOR_DIR}/clapack/include" )
    endif()

    if( NOT GSL_LIB_DIR AND IS_DIRECTORY $ENV{GSL_LIB_DIR}  )
        set( GSL_LIB_DIR $ENV{GSL_LIB_DIR} )
        set( GSL_INC_DIR $ENV{GSL_INC_DIR} )
    endif()
    if( NOT GSL_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/gsl/lib )
        set( GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib" )
        set( GSL_INC_DIR "${VENDOR_DIR}/gsl/include" )
    endif()

    set_package_properties( MPI PROPERTIES
       URL "http://www.open-mpi.org/"
       DESCRIPTION "A High Performance Message Passing Library"
       TYPE RECOMMENDED
       PURPOSE "If not available, all Draco components will be built as scalar applications."
       )
    set_package_properties( BLAS PROPERTIES
       DESCRIPTION "Basic Linear Algebra Subprograms"
       TYPE OPTIONAL
       PURPOSE "Required for bulding the lapack_wrap component." 
       )
    set_package_properties( LAPACK PROPERTIES
       DESCRIPTION "Linear Algebra PACKage"
       TYPE OPTIONAL
       PURPOSE "Required for bulding the lapack_wrap component." 
       )     
    set_package_properties( GSL PROPERTIES
       URL "http://www.gnu.org/s/gsl/"
       DESCRIPTION "GNU Scientific Library"
       TYPE REQUIRED
       PURPOSE "Required for bulding quadrature and rng components."
       )  

endmacro()

#------------------------------------------------------------------------------
# This macro should contain all the system libraries which are
# required to link the main objects.
#------------------------------------------------------------------------------
macro( setupVendorLibraries )

   message( "\nVendor Setup:\n")

  #
  # General settings
  # 
  setVendorVersionDefaults()

  # System specific settings
  if ( UNIX )
     setupMPILibrariesUnix()
#     if( "${SITE}" MATCHES "c[it]" )
        # # Provides BLAS and LAPACK
        # find_package( LIBSCI )
        # set_package_properties( LIBSCI PROPERTIES
        #    # URL "http://www.open-mpi.org/"
        #    DESCRIPTION "Cray's High Performance Scientify Library (LAPACK, BLAS, more)."
        #    TYPE RECOMMENDED
        #    PURPOSE 
        #    "Provides BLAS, LAPACK, BLACS, ScaLAPACK and SuperLU_DIST."
        #    )
#     else()
        setupLAPACKLibrariesUnix()
#     endif()
     setupVendorLibrariesUnix()
  elseif( WIN32 )
     setupVendorLibrariesWindows()
  else()
     message( FATAL_ERROR "
I don't know how to setup global (vendor) libraries for this platform.
WIN32=0; UNIX=0; CMAKE_SYSTEM=${CMAKE_SYSTEM}; 
CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}" )
  endif()

endmacro()

#----------------------------------------------------------------------#
# End vendor_libraries.cmake
#----------------------------------------------------------------------#
