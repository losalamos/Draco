#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-pgf90.cmake
# author Kelly Thompson 
# date   2011 June 7
# brief  Establish flags for Unix - PGI Fortran
# note   Copyright © 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "PGI" )

set( CMAKE_Fortran_FLAGS
  "-Mpreprocess"  ) 
set( CMAKE_Fortran_FLAGS_DEBUG 
  "-g -Mbounds -Mchkptr")
set( CMAKE_Fortran_FLAGS_RELEASE 
  "-O3")
set( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}" )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO 
  "${CMAKE_Fortran_FLAGS_DEBUG} -O3")

#------------------------------------------------------------------------------#
# End config/unix-pgf90.cmake
#------------------------------------------------------------------------------#
