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

set( DRACO_C4 "SCALAR" CACHE STRING 
   "Keyword for creating new libraries (SCALAR or MPI)."
   FORCE )

# clubimc/src/heterogeneous/tests/milagro/echo/milagro_rz_mg_Test_Host_Rep_accel_side.cc
# Host and Cell
option( HET_MESH_EVERY_CYCLE "See milagro_rz_mg_Test_Host_Rep_accel_side.cc" ON )
if( HET_MESH_EVERY_CYCLE )
   add_definitions( -DMESH_EVERY_CYCLE )
endif()

# Monitor requests and control messages from the host.
# Host and Cell
# clubimc/src/heterogeneous/accel_lib/DACS_Host.hh
# clubimc/src/heterogeneous/host/DACS_Accelerator.t.hh
# src/heterogeneous/host_accel/Particle_Comm_Buffer_dynamic.hh
# src/heterogeneous/ppe_lib/event_loop2.hh
option( HET_HOST_ACCEL_DACS " " ON )
if( HET_HOST_ACCEL_DACS )
   add_definitions( -DHOST_ACCEL_DACS )
endif()

# Host and cell
# clubimc/src/heterogeneous/host/Host_Particle_Rcvr.t.hh
# clubimc/src/heterogeneous/ppe_lib/Multiple_Particle_Reader.t.hh
option( HET_PPE_WRITE_BUFFER_DIRECT " " ON )
if( HET_PPE_WRITE_BUFFER_DIRECT )
   add_definitions( -DPPE_WRITE_BUFFER_DIRECT )
endif()

# Host and cell
# clubimc/src/heterogeneous/accel_lib/Accel_Particle_Rcvr.hh
# clubimc/src/heterogeneous/accel_lib/Accel_Particle_Rcvr.t.hh
# clubimc/src/heterogeneous/host/Host_Particle_Xmitter.t.hh
# clubimc/src/heterogeneous/ppe_lib/Multiple_Particle_Writer.i.hh
option( HET_PPE_READ_BUFFER_DIRECT  " " ON )
if( HET_PPE_READ_BUFFER_DIRECT )
   add_definitions( -DPPE_READ_BUFFER_DIRECT )
endif()



#--------------------------------------------------------------------------------

# Host only
# clubimc/src/heterogeneous/ppe_lib/event_loop2.t.hh
option( HET_ACCEL_RECV_IPROBE       " " OFF )
if( HET_ACCEL_RECV_IPROBE )
   add_definitions( -DACCEL_RECV_IPROBE )
endif()

if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )
   # Host and Cell
   # clubimc/src/heterogeneous/ppe_lib/event_loop2.t.hh
   option( HET_ACCEL_RECV_NONBLOCKING  " " ON )
   if( HET_ACCEL_RECV_NONBLOCKING )
      add_definitions( -DACCEL_RECV_NONBLOCKING )
   else()
      add_definitions( -DACCEL_RECV_BLOCKING )
   endif()

   # clubimc/src/heterogeneous/ppe_lib/event_loop2.t.hh
   # Host and Cell
   option( HET_ACCEL_SEND_NONBLOCKING  " " ON )
   if( HET_ACCEL_SEND_NONBLOCKING )
      add_definitions( -DACCEL_SEND_NONBLOCKING )
   else()
      add_definitions( -DACCEL_SEND_BLOCKING )
   endif()
else()
   # Host and Cell
   # clubimc/src/heterogeneous/ppe_lib/event_loop2.t.hh
   option( HET_HOST_RECV_NONBLOCKING  " " ON )
   if( HET_HOST_RECV_NONBLOCKING )
      add_definitions( -DHOST_RECV_NONBLOCKING )
   else()
      add_definitions( -DHOST_RECV_BLOCKING )
   endif()

   # clubimc/src/heterogeneous/ppe_lib/event_loop2.t.hh
   # Host and Cell
   option( HET_HOST_SEND_NONBLOCKING  " " ON )
   if( HET_HOST_SEND_NONBLOCKING )
      add_definitions( -DHOST_SEND_NONBLOCKING )
   else()
      add_definitions( -DHOST_SEND_BLOCKING )
   endif()
endif()


# Cell and SPE
# clubimc/src/heterogeneous/ppe_apps/accel_side_rz_mg/run_time_step.cc
# clubimc/src/heterogeneous/ppe_apps/accel_side_rz_mg/run_time_step.hh
# clubimc/src/heterogeneous/ppe_apps/accel_side_rz_mg/spe/run_particle_transporter.cc
# clubimc/src/heterogeneous/ppe_apps/accel_side_rz_mg/spe/run_particle_transporter_w_RW.cc
if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )
   option( HET_SHORT_SPE_TALLIES       " " ON )
   if( HET_SHORT_SPE_TALLIES )
      add_definitions( -DSHORT_SPE_TALLIES )   # this is for cell_cpp_flags
      # add_definitions( -DSHORT_SPE_TALLIES ) # must also be done for spe_cppflags.
   endif()
endif()

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )
   # On for x86 builds of ds++

   # Host only
   option( HET_HOST_DIRECT_DACS_INPUT " " ON )
   if( HET_HOST_DIRECT_DACS_INPUT )
      add_definitions( -DHOST_DIRECT_DACS_INPUT )
   endif()
   
   # Host only
   option( HET_HOST_DIRECT_DACS_OUTPUT " " OFF )
   if( HET_HOST_DIRECT_DACS_OUTPUT )
      add_definitions( -DHOST_DIRECT_DACS_OUTPUT )
   endif()
endif()


# Change some compiler flags for the roadrunner code
if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )
if( "${CMAKE_GENERATOR}" MATCHES "Makefiles" )
   message( "NOTICE: We are modifying the default g++ compile flags for roadrunner (clubimc/pkg_config/platform_customization_roadrunner_ppe.cmake)")

# -pthread -O0 -Wnon-virtual-dtor -Wreturn-type -pedantic
# -Wno-long-long -finline-functions -gdwarf-2 -DHOST_ACCEL_DACS
# -DMESH_EVERY_CYCLE -DPPE_WRITE_BUFFER_DIRECT
# -DPPE_READ_BUFFER_DIRECT -DHOST_SEND_NONBLOCKING
# -DHOST_RECV_NONBLOCKING -DHOST_DIRECT_DACS_INPUT -DOMPI_SKIP_MPICXX
# -DADDRESSING_64 -DDBC=7 

  set( DRACO_C_FLAGS                "-m64 -pthread -finline-functions -DADDRESSING_64" )
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