#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-xlf.cmake
# author Gabriel Rockefeller
# date   2012 Nov 1
# brief  Establish flags for Unix - IBM XL Fortran
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "XL" )

if( NOT Fortran_FLAGS_INITIALIZED )
   set( Fortran_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )
   set( CMAKE_Fortran_COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION} CACHE
        STRING "Fortran compiler version string" FORCE )
   mark_as_advanced( CMAKE_Fortran_COMPILER_VERSION )

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

#
# Toggle compiler flags for optional features
#
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "" )

# ----------------------------------------------------------------------------
# Helper macro to fix the value of
# CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES when compiling on SQ with
# xlc/xlf.  For some reason CMake is adding xlomp_ser to the link line
# when a C++ main is built by linking to F90 libraries.  This library
# should not be on For some reason CMake is adding xlomp_ser to the
# link line when a C++ main the link line and actually causes the link
# to fail on SQ.
# ----------------------------------------------------------------------------
macro( remove_lib_from_link lib_for_removal )

   set( tmp "" )
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

endmacro()

# When building on SQ with XLC/XLF90, ensure that libxlomp_ser is not
# on the link line
#remove_lib_from_link( "xlomp_ser" )

#------------------------------------------------------------------------------#
# End config/unix-xlf.cmake
#------------------------------------------------------------------------------#
