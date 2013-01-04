#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-xlf.cmake
# author Gabriel Rockefeller
# date   2012 Nov 1
# brief  Establish flags for Unix - IBM XL Fortran
# note   Copyright (C) 2012 Los Alamos National Security, LLC.
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

set( CMAKE_Fortran_FLAGS    "-qlanglvl=2003std -qinfo=all -qflag=i:w -qarch=auto" )
if( USE_OPENMP )
   set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp=omp" )
endif()

set( CMAKE_Fortran_FLAGS_DEBUG          "-g -O0 -qcheck" )
SET( CMAKE_Fortran_FLAGS_RELEASE        "-O3 -qhot=novector -qsimd=auto -qstrict=nans:operationprecision" )
SET( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE}" )
SET( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O3 -qhot=novector -qsimd=auto -qstrict=nans:operationprecision" )

#------------------------------------------------------------------------------#
# End config/unix-xlf.cmake
#------------------------------------------------------------------------------#
