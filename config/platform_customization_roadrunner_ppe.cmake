#-----------------------------*-cmake-*----------------------------------------#
# file   clubimc/pkg_config
# author Kelly Thompson <kgt@lanl.gov>
# date   2011 May 10
# brief  Custom build command for cross compiling Cell (PPE) code on roadrunner
# note   Copyright (C) 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# On RR default the Het code to 'ON'
option( ENABLE_HET_CODE "Build code for PPEs." ON )
add_feature_info( ENABLE_HET_CODE ENABLE_HET_CODE 
   "Build software for heterogeneous architectures (e.g. roadrunner).
   Required for building the device component.")

# Must use static libraries.
set( DRACO_LIBRARY_TYPE "STATIC" CACHE STRING 
   "Keyword for creating new libraries (STATIC or SHARED)."
   FORCE )

# Don't use MPI on the PPC or SPE
if( "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )
   set( DRACO_C4 "SCALAR" CACHE STRING 
      "Keyword for creating new libraries (SCALAR or MPI)."
      FORCE )
   unset( DRACO_SCALAR )
endif()

#----------------------------------------------------------------------
# Check requirements
#----------------------------------------------------------------------

# We have never been able to get the DaCS libraries to work correctly
# with pgCC.  On 9-13-2011, KT spent about 1 day trying to make this
# work to no avail.  The device tests will run as far as the
# DACS_Device destructor where the code reports
#
# ~DACS_Device reports (release): DACS_ERR_NOT_INITIALIZED
#
# and then the process hangs -- presumably due a corrupt DaCS state.
# 
if( "${ENABLE_HET_CODE}" AND "${CMAKE_CXX_COMPILER}" MATCHES "pgCC" )
   message( FATAL_ERROR "The DaCS libraries (cellsdk/3.1) are not "
      "compatible with the PGI C compiler." )
endif()


#----------------------------------------------------------------------
# Settings for the PPC and x86 (but not SPE)
#----------------------------------------------------------------------

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "spu-g[+][+]" )

   option( HET_MESH_EVERY_CYCLE 
      "See milagro_rz_mg_Test_Host_Rep_accel_side.cc" ON )
#   add_feature_info( HET_MESH_EVERY_CYCLE HET_MESH_EVERY_CYCLE 
#   "Possibly deprecated?")
   if( HET_MESH_EVERY_CYCLE )
      add_definitions( -DMESH_EVERY_CYCLE )
   endif()
   
   option( HET_HOST_ACCEL_DACS " " ON )
   add_feature_info( HET_HOST_ACCEL_DACS HET_HOST_ACCEL_DACS
      "No sure of purpose.")
   if( HET_HOST_ACCEL_DACS )
      add_definitions( -DHOST_ACCEL_DACS )
   endif()
   
   option( HET_PPE_WRITE_BUFFER_DIRECT " " ON )
   add_feature_info( HET_PPE_WRITE_BUFFER_DIRECT
   HET_PPE_WRITE_BUFFER_DIRECT 
   "There are 2 implementaitons of drain_buffer.  When enabled, the PPE sends particles via buffers as used on the accelerator instead of sending packed particles (Particle::pack()).")
   if( HET_PPE_WRITE_BUFFER_DIRECT )
      add_definitions( -DPPE_WRITE_BUFFER_DIRECT )
   endif()
   
   option( HET_PPE_READ_BUFFER_DIRECT  " " ON )
   add_feature_info( HET_PPE_READ_BUFFER_DIRECT
   HET_PPE_READ_BUFFER_DIRECT 
   "There are 2 implementaitons of drain_buffer.  When enabled, the PPE sends particles via buffers as used on the accelerator instead of sending packed particles (Particle::pack()).")
   if( HET_PPE_READ_BUFFER_DIRECT )
      add_definitions( -DPPE_READ_BUFFER_DIRECT )
   endif()

   option( HET_ACCEL_RECV_IPROBE       " " OFF )
   add_feature_info( HET_ACCEL_RECV_IPROBE HET_ACCEL_RECV_IPROBE
   "Unknown.  Possibly deprecated, see ppe_lib/event_loop2.t.hh")
   if( HET_ACCEL_RECV_IPROBE )
      add_definitions( -DACCEL_RECV_IPROBE )
   endif()

   option( HET_HOST_ACCEL_DACS_GROUP " " ON )
   add_feature_info( HET_HOST_ACCEL_DACS_GROUP HET_HOST_ACCEL_DACS_GROUP "Unknown.")
   if( HET_HOST_ACCEL_DACS_GROUP )
      add_definitions( -DHOST_ACCEL_DACS_GROUP )
   endif()

endif()

#----------------------------------------------------------------------
# Settings for the PPC (not x86 or SPE)
#----------------------------------------------------------------------

if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )

   option( HET_ACCEL_RECV_NONBLOCKING  " " ON )
   add_feature_info( HET_ACCEL_RECV_NONBLOCKING
      HET_ACCEL_RECV_NONBLOCKING 
      "Use non-blocking receives in event_loop2.t.hh")
   if( HET_ACCEL_RECV_NONBLOCKING )
      add_definitions( -DACCEL_RECV_NONBLOCKING )
   else()
      add_definitions( -DACCEL_RECV_BLOCKING )
   endif()

   option( HET_ACCEL_SEND_NONBLOCKING  " " ON )
   add_feature_info( HET_ACCEL_SEND_NONBLOCKING
      HET_ACCEL_SEND_NONBLOCKING 
      "Use non-blocking sends in event_loop2.t.hh")
   if( HET_ACCEL_SEND_NONBLOCKING )
      add_definitions( -DACCEL_SEND_NONBLOCKING )
   else()
      add_definitions( -DACCEL_SEND_BLOCKING )
   endif()

   option( HET_CHECK_PARTICLES_PPE_RECV " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_PPE_RECV
      HET_CHECK_PARTICLES_PPE_RECV 
      "Perform extra validity checks on particle objects on the PPE.")
   if( HET_CHECK_PARTICLES_PPE_RECV )
      add_definitions( -DCHECK_PARTICLES_PPE_RECV )
   endif()

   option( HET_CHECK_PARTICLES_PPE_SEND " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_PPE_SEND
      HET_CHECK_PARTICLES_PPE_SEND 
      "Perform extra validity checks on particle objects on the PPE.")
   if( HET_CHECK_PARTICLES_PPE_SEND )
      add_definitions( -DCHECK_PARTICLES_PPE_SEND )
   endif()

endif()

#----------------------------------------------------------------------
# Settings for the x86 (not PPE or SPE)
#----------------------------------------------------------------------

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )

   option( HET_HOST_RECV_NONBLOCKING  " " ON )
   add_feature_info( HET_HOST_RECV_NONBLOCKING
      HET_HOST_RECV_NONBLOCKING
      "Unknown.")
   if( HET_HOST_RECV_NONBLOCKING )
      add_definitions( -DHOST_RECV_NONBLOCKING )
   else()
      add_definitions( -DHOST_RECV_BLOCKING )
   endif()

   option( HET_HOST_SEND_NONBLOCKING  " " ON )
   add_feature_info( HET_HOST_SEND_NONBLOCKING
      HET_HOST_SEND_NONBLOCKING
      "Unknown.")
   if( HET_HOST_SEND_NONBLOCKING )
      add_definitions( -DHOST_SEND_NONBLOCKING )
   else()
      add_definitions( -DHOST_SEND_BLOCKING )
   endif()

   option( HET_HOST_DIRECT_DACS_INPUT " " ON )
   add_feature_info( HET_HOST_DIRECT_DACS_INPUT
      HET_HOST_DIRECT_DACS_INPUT
      "Unknown.")
   if( HET_HOST_DIRECT_DACS_INPUT )
      add_definitions( -DHOST_DIRECT_DACS_INPUT )
   endif()
   
   option( HET_HOST_DIRECT_DACS_OUTPUT " " OFF )
   add_feature_info( HET_HOST_DIRECT_DACS_OUTPUT
      HET_HOST_DIRECT_DACS_OUTPUT
      "Unknown.")
   if( HET_HOST_DIRECT_DACS_OUTPUT )
      add_definitions( -DHOST_DIRECT_DACS_OUTPUT )
   endif()

   option( HET_CHECK_PARTICLES_HOST_DISPATCH " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_HOST_DISPATCH
      HET_CHECK_PARTICLES_HOST_DISPATCH
      "Perform extra validity checks on paticle objects.")
   if( HET_CHECK_PARTICLES_HOST_DISPATCH )
      add_definitions( -DCHECK_PARTICLES_HOST_DISPATCH )
   endif()

   option( HET_CHECK_PARTICLES_HOST_BANK_RECV " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_HOST_BANK_RECV
      HET_CHECK_PARTICLES_HOST_BANK_RECV
      "Perform extra validity checks on paticle objects.")
   if( HET_CHECK_PARTICLES_HOST_BANK_RECV )
      add_definitions( -DCHECK_PARTICLES_HOST_BANK_RECV )
   endif()

   option( HET_CHECK_PARTICLES_HOST_STREAM_SOURCE " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_HOST_STREAM_SOURCE
      HET_CHECK_PARTICLES_HOST_STREAM_SOURCE
      "Perform extra validity checks on paticle objects.")
   if( HET_CHECK_PARTICLES_HOST_STREAM_SOURCE )
      add_definitions( -DCHECK_PARTICLES_HOST_STREAM_SOURCE )
   endif()

   option( HET_CHECK_PARTICLES_HOST_VOL_SOURCE " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_HOST_VOL_SOURCE
      HET_CHECK_PARTICLES_HOST_VOL_SOURCE
      "Perform extra validity checks on paticle objects.")
   if( HET_CHECK_PARTICLES_HOST_VOL_SOURCE )
      add_definitions( -DCHECK_PARTICLES_HOST_VOL_SOURCE )
   endif()

   option( HET_CHECK_PARTICLES_HOST_CEN_SOURCE " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_HOST_CEN_SOURCE
      HET_CHECK_PARTICLES_HOST_CEN_SOURCE
      "Perform extra validity checks on paticle objects.")
   if( HET_CHECK_PARTICLES_HOST_CEN_SOURCE )
      add_definitions( -DCHECK_PARTICLES_HOST_CEN_SOURCE )
   endif()

endif()

#----------------------------------------------------------------------
# Settings for the PPE and SPE (but not for x86)
#----------------------------------------------------------------------

if( "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )
   option( HET_SHORT_SPE_TALLIES       " " ON )
   add_feature_info( HET_SHORT_SPE_TALLIES
      HET_SHORT_SPE_TALLIES
      "If enabled use rtt_comm::Private_Short_Tally_Receiver<Tally_T, NSPES> in place of rtt_comm::Private_Tally_Receiver<Tally_T, NSPES> .")
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
   add_feature_info( HET_FLAT_AMR_RDR
      HET_FLAT_AMR_RDR
      "Toggle between the normal and a flattened AMR mesh reader.")
   if( HET_FLAT_AMR_RDR )
      add_definitions( -DFLAT_AMR_RDR )
   endif()

   option( HET_FLAT_OP_RDR " " ON )
   add_feature_info( HET_FLAT_OP_RDR
      HET_FLAT_OP_RDR
      "Toggle between the normal and a flattened opacity reader.")
   if( HET_FLAT_OP_RDR )
      add_definitions( -DFLAT_OP_RDR )
   endif()

   option( HET_CHECK_PARTICLES_SPE_RECV " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_SPE_RECV
      HET_CHECK_PARTICLES_SPE_RECV
      "Perform extra validity checks on paticle objects.")
   if( HET_CHECK_PARTICLES_SPE_RECV )
      add_definitions( -DCHECK_PARTICLES_SPE_RECV )
   endif()

   option( HET_CHECK_PARTICLES_SPE_SEND " " OFF )
   add_feature_info( HET_CHECK_PARTICLES_SPE_SEND
      HET_CHECK_PARTICLES_SPE_SEND
      "Perform extra validity checks on paticle objects.")
   if( HET_CHECK_PARTICLES_SPE_SEND )
      add_definitions( -DCHECK_PARTICLES_SPE_SEND )
   endif()
endif()

##---------------------------------------------------------------------------##
## End platform_customization_roadrunner_ppe.cmake
##---------------------------------------------------------------------------##

