#-----------------------------*-cmake-*----------------------------------------#
# file   clubimc/pkg_config
# author Kelly Thompson <kgt@lanl.gov>
# date   2011 May 10
# brief  Custom build command for cross compiling Cell (PPE) code on roadrunner
# note   Copyright © 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#


# Must use static libraries.
set( DRACO_LIBRARY_TYPE "STATIC" CACHE STRING 
   "Keyword for creating new libraries (STATIC or SHARED)."
   FORCE )

if( "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )
   set( DRACO_C4 "SCALAR" CACHE STRING 
      "Keyword for creating new libraries (SCALAR or MPI)."
      FORCE )
   unset( DRACO_SCALAR )
endif()

#----------------------------------------------------------------------
# Settings for the PPC and x86 (but not SPE)
#----------------------------------------------------------------------

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "spu-g[+][+]" )

   option( HET_MESH_EVERY_CYCLE "See milagro_rz_mg_Test_Host_Rep_accel_side.cc" ON )
   if( HET_MESH_EVERY_CYCLE )
      add_definitions( -DMESH_EVERY_CYCLE )
   endif()
   
   option( HET_HOST_ACCEL_DACS " " ON )
   if( HET_HOST_ACCEL_DACS )
      add_definitions( -DHOST_ACCEL_DACS )
   endif()
   
   option( HET_PPE_WRITE_BUFFER_DIRECT " " ON )
   if( HET_PPE_WRITE_BUFFER_DIRECT )
      add_definitions( -DPPE_WRITE_BUFFER_DIRECT )
   endif()
   
   option( HET_PPE_READ_BUFFER_DIRECT  " " ON )
   if( HET_PPE_READ_BUFFER_DIRECT )
      add_definitions( -DPPE_READ_BUFFER_DIRECT )
   endif()

   option( HET_ACCEL_RECV_IPROBE       " " OFF )
   if( HET_ACCEL_RECV_IPROBE )
      add_definitions( -DACCEL_RECV_IPROBE )
   endif()

   option( HET_HOST_ACCEL_DACS_GROUP " " ON )
   if( HET_HOST_ACCEL_DACS_GROUP )
      add_definitions( -DHOST_ACCEL_DACS_GROUP )
   endif()

endif()

#----------------------------------------------------------------------
# Settings for the PPC (not x86 or SPE)
#----------------------------------------------------------------------

if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )

   option( HET_ACCEL_RECV_NONBLOCKING  " " ON )
   if( HET_ACCEL_RECV_NONBLOCKING )
      add_definitions( -DACCEL_RECV_NONBLOCKING )
   else()
      add_definitions( -DACCEL_RECV_BLOCKING )
   endif()

   option( HET_ACCEL_SEND_NONBLOCKING  " " ON )
   if( HET_ACCEL_SEND_NONBLOCKING )
      add_definitions( -DACCEL_SEND_NONBLOCKING )
   else()
      add_definitions( -DACCEL_SEND_BLOCKING )
   endif()

endif()

#----------------------------------------------------------------------
# Settings for the x86 (not PPE or SPE)
#----------------------------------------------------------------------

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )

   option( HET_HOST_RECV_NONBLOCKING  " " ON )
   if( HET_HOST_RECV_NONBLOCKING )
      add_definitions( -DHOST_RECV_NONBLOCKING )
   else()
      add_definitions( -DHOST_RECV_BLOCKING )
   endif()

   option( HET_HOST_SEND_NONBLOCKING  " " ON )
   if( HET_HOST_SEND_NONBLOCKING )
      add_definitions( -DHOST_SEND_NONBLOCKING )
   else()
      add_definitions( -DHOST_SEND_BLOCKING )
   endif()

   option( HET_HOST_DIRECT_DACS_INPUT " " ON )
   if( HET_HOST_DIRECT_DACS_INPUT )
      add_definitions( -DHOST_DIRECT_DACS_INPUT )
   endif()
   
   option( HET_HOST_DIRECT_DACS_OUTPUT " " OFF )
   if( HET_HOST_DIRECT_DACS_OUTPUT )
      add_definitions( -DHOST_DIRECT_DACS_OUTPUT )
   endif()

endif()

#----------------------------------------------------------------------
# Settings for the PPE and SPE (but not for x86)
#----------------------------------------------------------------------

if( "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )
   option( HET_SHORT_SPE_TALLIES       " " ON )
   if( HET_SHORT_SPE_TALLIES )
      add_definitions( -DSHORT_SPE_TALLIES )   # this is for cell_cpp_flags
      # add_definitions( -DSHORT_SPE_TALLIES ) # must also be done for spe_cppflags.
   endif()
endif()

#----------------------------------------------------------------------
# Settings for the SPE (not x86 or PPE)
#----------------------------------------------------------------------

if( "${CMAKE_CXX_COMPILER}" MATCHES "spu-g[+][+]" )
   option( HET_FLAT_AMR_RDR " " ON )
   if( HET_FLAT_AMR_RDR )
      add_definitions( -DFLAT_AMR_RDR )
   endif()

   option( HET_FLAT_OP_RDR " " ON )
   if( HET_FLAT_OP_RDR )
      add_definitions( -DFLAT_OP_RDR )
   endif()
endif()

#----------------------------------------------------------------------
# Use Tim Kelley's CXX flags for g++ instead of the default draco
# flags. 
# 
# \BUG: This section should be removed after testing with the default
# flags is complete/resolved.
#----------------------------------------------------------------------

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )
   if( "${CMAKE_GENERATOR}" MATCHES "Makefiles" )
      message( "NOTICE: We are modifying the default g++ compile flags for roadrunner (clubimc/pkg_config/platform_customization_roadrunner_ppe.cmake)")

      set( DRACO_C_FLAGS                "-m64 -pthread -finline-functions -DADDRESSING_64 -DCACHE_LINE_SIZE=128" )
      set( DRACO_C_FLAGS_DEBUG          "-gdwarf-2 -O0 -DDEBUG")
      set( DRACO_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
      set( DRACO_C_FLAGS_MINSIZEREL     "${DRACO_C_FLAGS_RELEASE}" )
      set( DRACO_C_FLAGS_RELWITHDEBINFO "${DRACO_C_FLAGS_DEBUG} -O3 -DDEBUG" )
      
      set( DRACO_CXX_FLAGS                "${DRACO_C_FLAGS}" )
      set( DRACO_CXX_FLAGS_DEBUG          "${DRACO_C_FLAGS_DEBUG} -pedantic -Wnon-virtual-dtor -Wreturn-type -Woverloaded-virtual -Wno-long-long")
      set( DRACO_CXX_FLAGS_RELEASE        "${DRACO_C_FLAGS_RELEASE}")
      set( DRACO_CXX_FLAGS_MINSIZEREL     "${DRACO_CXX_FLAGS_RELEASE}")
      set( DRACO_CXX_FLAGS_RELWITHDEBINFO "${DRACO_C_FLAGS_RELWITHDEBINFO}" )
      
      if( "${HET_CXX_FLAGS_INITIALIZED}no" STREQUAL "no" )
         set( HET_CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL 
            "using draco settings (het)." )
         set( CMAKE_C_FLAGS                "${DRACO_C_FLAGS}"                CACHE STRING "compiler flags" FORCE )
         set( CMAKE_C_FLAGS_DEBUG          "${DRACO_C_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
         set( CMAKE_C_FLAGS_RELEASE        "${DRACO_C_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
         set( CMAKE_C_FLAGS_MINSIZEREL     "${DRACO_C_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
         set( CMAKE_C_FLAGS_RELWITHDEBINFO "${DRACO_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )
         set( CMAKE_CXX_FLAGS                "${DRACO_CXX_FLAGS}"                CACHE STRING "compiler flags" FORCE )
         set( CMAKE_CXX_FLAGS_DEBUG          "${DRACO_CXX_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
         set( CMAKE_CXX_FLAGS_RELEASE        "${DRACO_CXX_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
         set( CMAKE_CXX_FLAGS_MINSIZEREL     "${DRACO_CXX_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
         set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${DRACO_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )
      endif()
      
   endif()
endif()


#----------------------------------------------------------------------
# Instructions for building the SPE libraries.
#
# This is done as part of the PPE build at 'make-time' not at
# configure-time
#----------------------------------------------------------------------
macro( ClubIMC_build_spe_libraries)
   if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )
      message( "
-- SPE libraries will be compiled:
     SPE binary directory : ${ClubIMC_BINARY_DIR}/clubimc_spe
     SPE install directory: ${CMAKE_INSTALL_PREFIX}/spe
")
      include(ExternalProject)
      ExternalProject_Add( ClubIMC_SPE
         # Root dir for build
         PREFIX            ${ClubIMC_BINARY_DIR}/clubimc_spe
         SOURCE_DIR        ${ClubIMC_SOURCE_DIR}
         CONFIGURE_COMMAND cmake
                        -DCMAKE_TOOLCHAIN_FILE=${DRACO_DIR}/config/Toolchain-roadrunner-spu.cmake 
                        -DDRACO_DIR=${DRACO_DIR} 
                        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/spe
                        ${ClubIMC_SOURCE_DIR}
         CMAKE_GENERATOR   ${CMAKE_GENERATOR}
         BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM}
         INSTALL_COMMAND   ${CMAKE_MAKE_PROGRAM} install
         )
      set( ClubIMC_SPE_DIR ${CMAKE_INSTALL_PREFIX}/spe )
      set( ClubIMC_SPE_rz_mg_11_spe_LIBRARY ${ClubIMC_SPE_DIR}/lib/librtt_rz_mg_11_spe.a )
      set( ClubIMC_SPE_rz_mg_11_spe_rw_LIBRARY ${ClubIMC_SPE_DIR}/lib/librtt_rz_mg_11_spe_rw.a )
   endif()
endmacro()

