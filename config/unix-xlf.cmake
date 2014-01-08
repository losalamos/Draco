#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-xlf.cmake
# author Gabriel Rockefeller
# date   2012 Nov 1
# brief  Establish flags for Unix - IBM XL Fortran
# note   Copyright (C) 2012-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "XL" )

# Notable xlf options:

# -qlanglvl=2003std
#                 Accept ISO Fortran 2003 standard language features.
# -qinfo=all      Enable all warnings (and informational messages).
# -qflag=i:w      Send informational-level messages to a listing file,
#                 if one is requested, but only send warning-level or
#                 more severe messages (i.e., errors) to the terminal.
# -qcheck         Enable array element, array section, and character
#                 substring checks.
# -O3             -O2 plus memory- and compile-time-intensive
#                 operations that can alter the semantics of programs.
# -qstrict=nans:operationprecision
#                 Disable optimizations at -O3 and above that may
#                 produce incorrect results in the presence of NaNs,
#                 or that produce approximate results for individual
#                 floating-point operations.
# -qhot=novector  Enable high-order transformations during
#                 optimization.  LLNL recommends novector to disable
#                 gathering math intrinsics into separate vector math
#                 library calls (because it's typically better to let
#                 those instructions intermix with other
#                 floating-point operations, when using SIMD
#                 instructions).
# -qsimd=auto     Enable automatic generation of SIMD instructions, to
#                 take advantage of BG/Q-specific Quad Processing
#                 eXtension (QPX) units.
# -qsmp=omp       Enable parallelization using OpenMP pragmas.

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_Fortran_FLAGS                "-qlanglvl=2003std -qinfo=all -qflag=i:w -qarch=auto" )
   set( CMAKE_Fortran_FLAGS_DEBUG          "-g -O0 -qnosmp -qcheck" )
   set( CMAKE_Fortran_FLAGS_RELEASE        "-O3 -qhot=novector -qsimd=auto -qstrict=nans:operationprecision" )
   set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE}" )
   set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O3 -qhot=novector -qsimd=auto -qstrict=nans:operationprecision" )

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

toggle_compiler_flag( USE_OPENMP "-qsmp=omp" "Fortran" "")

# ----------------------------------------------------------------------------
# Helper macro to fix the value of
# CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES when compiling on SQ with
# xlc/xlf.  For some reason CMake is adding xlomp_ser to the link line
# when a C++ main is built by linking to F90 libraries.  This library
# should not be on For some reason CMake is adding xlomp_ser to the
# link line when a C++ main the link line and actually causes the link
# to fail on SQ.
# ----------------------------------------------------------------------------
function( remove_lib_from_link lib_for_removal )

   foreach( lib ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )
      if( ${lib} MATCHES ${lib_for_removal} )
         message("Removing ${lib_for_removal} from CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES")
         # do not add ${lib_for_removal} to the list ${tmp}.
      else()
         # Add a libraries except for the above match to the new list.
         list( APPEND tmp ${lib} )
      endif()
   endforeach()

   # Ensure that the value is saved to CMakeCache.txt
   set( CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "${tmp}" CACHE STRING
      "Extra Fortran libraries needed when linking to Fortral library from C++ target"
      FORCE )

endfunction()

# When building on SQ with XLC/XLF90, ensure that libxlomp_ser is not
# on the link line
remove_lib_from_link( "xlomp_ser" )

#------------------------------------------------------------------------------#
# End config/unix-xlf.cmake
#------------------------------------------------------------------------------#
