#-----------------------------*-cmake-*----------------------------------------#
# file   src/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building src level Makefile.
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( draco_src_dir )

# Platform Checks 
# Is HOST_NAME_MAX defined?  Is WinSock2.h available?  Is gethostname() available?
message("
Platform Checks...
")
include( platform_checks )
set_draco_uname()
query_have_gethostname()
query_have_maxpathlen()
query_have_sys_headers() # sets HAVE_UNISTD_H, etc.
query_have_restrict_keyword()
if( DRACO_ENABLE_CXX11 )
   set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
      "${Draco_SOURCE_DIR}/config/cxx11")
   include(CheckCXX11Features)
   # Results are saved in ds++/config.h
endif()

# Provide helper functions used by component CMakeLists.txt files
include( component_macros )

# Helper macro that only tries to configure a listed package if it
# exists in the source tree.  For example, when building milagro, the
# get_draco script only checks out the parts of Draco that are
# required by milagro.
macro( add_dir_if_exists package )
   message( "   ${package}" )
   if( EXISTS ${PROJECT_SOURCE_DIR}/${package} )
      add_subdirectory( ${package} )
   endif()
endmacro()

# Level 1
message(" ")
message( STATUS "Configuring Level 1 packages --" )
add_dir_if_exists( ds++ )

# Level 2
message(" ")
message( STATUS "Configuring Level 2 packages --" )
add_dir_if_exists( c4 )           # needs ds++
add_dir_if_exists( cdi )          # needs ds++
add_dir_if_exists( device )       # needs ds++
add_dir_if_exists( fpe_trap )     # needs ds++
if (NOT APPLE)
  add_dir_if_exists( lapack_wrap )  # needs ds++
endif()
add_dir_if_exists( linear )       # needs ds++
add_dir_if_exists( mesh_element ) # needs ds++
add_dir_if_exists( ode )          # needs ds++
add_dir_if_exists( rng )          # needs ds++
add_dir_if_exists( traits )       # needs ds++
add_dir_if_exists( units )        # needs ds++
if( _LANGUAGES_ MATCHES Fortran )
   if( NOT ${CMAKE_GENERATOR} MATCHES Xcode)
      add_dir_if_exists( FortranChecks )
   endif()
endif()

if( GRACE_FOUND )
  add_dir_if_exists( plot2D )
endif()

# Level 3
message(" ")
message( STATUS "Configuring Level 3 packages --" )
add_dir_if_exists( diagnostics )  # needs c4
add_dir_if_exists( fit )          # needs linear
add_dir_if_exists( meshReaders )  # needs c4
add_dir_if_exists( min )          # needs linear
add_dir_if_exists( norms )        # needs c4
add_dir_if_exists( parser )       # needs units, c4
add_dir_if_exists( roots )        # needs linear
add_dir_if_exists( timestep )     # needs c4, ds++
add_dir_if_exists( viz )          # needs traits, ds++
add_dir_if_exists( cdi_ipcress )

# Level 4
message(" ")
message( STATUS "Configuring Level 4 packages --" )
add_dir_if_exists( special_functions ) # needs roots, units, ode
add_dir_if_exists( cdi_analytic )      # needs parser, ode, cdi
add_dir_if_exists( cdi_eospac )        # needs parser, ode, cdi
add_dir_if_exists( RTT_Format_Reader)  # needs meshReaders
add_dir_if_exists( memory )            # needs diagnostics

# Level 5
if( NOT WIN32 )
  # this is broken on win32 due to a problem with Abstract_Class_Parser
  message(" ")
  message( STATUS "Configuring Level 5 packages --" )
  add_dir_if_exists( quadrature )       # needs mesh_element, parser, special_functions
endif()

# Summary
message(" ")

add_feature_info( BUILD_AUTODOC BUILD_AUTODOC
   "Turn on to enable generation of doxygen documentation (make autodoc).")

feature_summary( 
   WHAT PACKAGES_FOUND 
   INCLUDE_QUIET_PACKAGES
   DESCRIPTION "Packages found:")

feature_summary( 
   WHAT PACKAGES_NOT_FOUND 
   INCLUDE_QUIET_PACKAGES
   FATAL_ON_MISSING_REQUIRED_PACKAGES
   DESCRIPTION "Packages missing:")

feature_summary( 
   WHAT ENABLED_FEATURES 
   INCLUDE_QUIET_PACKAGES
   DESCRIPTION "Enabled features:")

feature_summary( 
   WHAT DISABLED_FEATURES 
   INCLUDE_QUIET_PACKAGES
   DESCRIPTION "Disabled features:")

message("Draco build summary:

Version     : Draco ${DRACO_VERSION_FULL}
Build type  : ${CMAKE_BUILD_TYPE}
C4 Model    : ${DRACO_C4}
Prefix dir  : ${CMAKE_INSTALL_PREFIX}
Source dir  : ${Draco_SOURCE_DIR}
Build dir   : ${Draco_BINARY_DIR}

CXX Compiler: ${CMAKE_CXX_COMPILER}
CXX FLAGS   : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}
C FLAGS     : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
if( _LANGUAGES_ MATCHES Fortran) 
   message("Fortran     : ${CMAKE_Fortran_COMPILER}")
   message("FC FLAGS    : ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
endif()
message(
"mpirun cmd  : ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} N ${MPIEXEC_POSTFLAGS_STRING}")
message("Library Type: ${DRACO_LIBRARY_TYPE}

ENABLE_RNG_NR    : ${ENABLE_RNG_NR}
DRACO_DBC_LEVEL  : ${DRACO_DBC_LEVEL}")
message("DRACO_DIAGNOSTICS: ${DRACO_DIAGNOSTICS}
DRACO_TIMING     : ${DRACO_TIMING}")

# Push some variables up one level
set( Draco_EXPORT_TARGET_PROPERTIES
   "${Draco_EXPORT_TARGET_PROPERTIES}" PARENT_SCOPE)

##---------------------------------------------------------------------------##
## Debug build system state:
##---------------------------------------------------------------------------##
option( DBS_PRINT_DEBUG_INFO "Print extra build system data" OFF )
if( DBS_PRINT_DEBUG_INFO )
   include( debug_macros )
   # Debug variables
   message( "CMake variables and values:" )
   echo_all_cmake_variable_values()

   # Debug target properties
   message( "CMake targets and properties:" )
   set( dbs_targets
      # Lib_dsxx
      # Lib_dsxx_test
      # Lib_c4
      # Lib_c4_test
      # Lib_cdi
      # Lib_cdi_test
      # Lib_fpe_trap
      # Lib_linear
      # Lib_linear_test
      # Lib_mesh_element
      # Lib_mesh_element_test
      # Lib_ode
      # Lib_ode_test
      # Lib_rng
      # Lib_rng_test
      # Lib_units
      # Lib_units_test
      # Lib_FC_Derived_Type
      # Lib_FortranChecks
      # Lib_plot2D
      # Lib_diagnostics
      # Lib_diagnostics_test
      # Lib_fit
      # Lib_fit_test
      # Lib_meshReaders
      # Lib_meshReaders_test
      # Lib_min
      # Lib_min_test
      # Lib_norms
      # Lib_norms_test
      # Lib_parser
      # Lib_parser_test
      # Lib_roots
      # Lib_roots_test
      # Lib_timestep
      # Lib_timestep_test
      # Lib_viz
      # Lib_viz_test
      # Lib_cdi_ipcress
      # Lib_cdi_ipcress_test
      # Lib_special_functions
      # Lib_cdi_analytic
      # Lib_cdi_analytic_test
      # Lib_RTT_Format_Reader
      # Lib_quadrature
      # Lib_quadrature_test
      # Ut_dsxx_tstArray_exe
      # Ut_dsxx_tstAssert_exe
      # Ut_dsxx_tstBounds_exe
      # Ut_dsxx_tstCheck_Strings_exe
      # Ut_dsxx_tstDBC_Array_exe
      # Ut_dsxx_tstData_Table_exe
      # Ut_dsxx_tstEndian_exe
      # Ut_dsxx_tstField_Traits_exe
      # Ut_dsxx_tstFile_Streams_exe
      # Ut_dsxx_tstHomogeneous_New_exe
      # Ut_dsxx_tstIndex_Converter_exe
      # Ut_dsxx_tstIndex_Counter_exe
      # Ut_dsxx_tstIndex_Set_exe
      # Ut_dsxx_tstMat1RA_exe
      # Ut_dsxx_tstMat2RA_exe
      # Ut_dsxx_tstMat3RA_exe
      # Ut_dsxx_tstMat4RA_exe
      # Ut_dsxx_tstMat5RA_exe
      # Ut_dsxx_tstMat_exe
      # Ut_dsxx_tstPacking_Utils_exe
      # Ut_dsxx_tstPath_exe
      # Ut_dsxx_tstRCF_exe
      # Ut_dsxx_tstRange_finder_exe
      # Ut_dsxx_tstSP_exe
      # Ut_dsxx_tstSafe_Divide_exe
      # Ut_dsxx_tstScalarUnitTest_exe
      # Ut_dsxx_tstSlice_exe
      # Ut_dsxx_tstSoft_Equiv_exe
      # Ut_dsxx_tstSortPermutation_exe
      # Ut_dsxx_tstVector_Lite_exe
      # Ut_dsxx_tstabs_exe
      # Ut_dsxx_tstconj_exe
      # Ut_dsxx_tstcube_exe
      # Ut_dsxx_tstdbc_exe
      # Ut_dsxx_tstisFinite_exe
      # Ut_dsxx_tstpythag_exe
      # Ut_dsxx_tstsign_exe
      # Ut_dsxx_tstsquare_exe
      # Ut_dsxx_tstto_string_exe
      # Ut_dsxx_tstvalue_exe
      # Ut_c4_phw_exe
      # Ut_c4_tstAbort_exe
      # Ut_c4_tstApplicationUnitTest_exe
      # Ut_c4_tstBroadcast_exe
      # Ut_c4_tstC4_Req_exe
      # Ut_c4_tstComm_Dup_exe
      # Ut_c4_tstCompare_exe
      # Ut_c4_tstGatherScatter_exe
      # Ut_c4_tstInvert_Comm_Map_exe
      # Ut_c4_tstOMP_exe
      # Ut_c4_tstParallelUnitTestFailMode_exe
      # Ut_c4_tstParallelUnitTest_exe
      # Ut_c4_tstPingPong_exe
      # Ut_c4_tstProcessor_Group_exe
      # Ut_c4_tstReduction_exe
      # Ut_c4_tstScalar_exe
      # Ut_c4_tstSend_Receive_exe
      # Ut_c4_tstSwap_exe
      # Ut_c4_tstTermination_Detector_exe
      # Ut_c4_tstTime_exe
      # Ut_c4_tstglobal_containers_exe
      # Ut_cdi_tCDI_exe
      # Ut_cdi_tDummyEoS_exe
      # Ut_cdi_tDummyOpacity_exe
      # Ut_fpe_trap_do_exception_exe
      # Ut_lapack_wrap_tstBlas_Level_1_exe
      # Ut_linear_tstbtridag_exe
      # Ut_linear_tstfnorm_exe
      # Ut_linear_tstgaussj_exe
      # Ut_linear_tstludcmp_exe
      # Ut_linear_tstqr_unpack_exe
      # Ut_linear_tstqrdcmp_exe
      # Ut_linear_tstqrupdt_exe
      # Ut_linear_tstrotate_exe
      # Ut_linear_tstrsolv_exe
      # Ut_linear_tstsvbksb_exe
      # Ut_linear_tstsvdcmp_exe
      # Ut_linear_tsttqli_exe
      # Ut_linear_tsttred2_exe
      # Ut_mesh_element_TestElementDefinition_exe
      # Ut_ode_tstFunction_Traits_exe
      # Ut_ode_tstquad_exe
      # Ut_ode_tstrkqs_exe
      # Ut_rng_tstRnd_Control_Inline_exe
      # Ut_rng_tstSubrandom_Sequence_exe
      # Ut_traits_tstViz_Traits_exe
      # Ut_units_tstFundUnit_exe
      # Ut_units_tstPhysicalConstants_exe
      # Ut_units_tstUnitSystemEnums_exe
      # Ut_units_tstUnitSystemType_exe
      # Ut_units_tstUnitSystem_exe
      # Ut_FortranChecks_cppmain_exe
      # Ut_FortranChecks_ftest_derived_types_exe
      # Ut_FortranChecks_ftest_f90main_exe
      # Ut_plot2D_install_inputs_1
      # Ut_plot2D_tstPlot2D_exe
      # Ut_diagnostics_tstDiagnostics_exe
      # Ut_diagnostics_tstTiming_exe
      # Ut_fit_tstsvdfit_exe
      # Ut_meshReaders_TestHexMeshReader_exe
      # Ut_meshReaders_install_inputs_1
      # Ut_min_install_inputs_1
      # Ut_min_tstbrent_exe
      # Ut_min_tstmnbrak_exe
      # Ut_min_tstmrqmin_exe
      # Ut_min_tstpowell_exe
      # Ut_norms_tstL2norm_exe
      # Ut_norms_tst_Norms_exe
      # Ut_parser_driver4tstConsole_Token_Stream_exe
      # Ut_parser_install_inputs_1
      # Ut_parser_tstAbstract_Class_Parser_exe
      # Ut_parser_tstClass_Parser_exe
      # Ut_parser_tstConsole_Token_Stream
      # Ut_parser_tstExpression_exe
      # Ut_parser_tstFile_Token_Stream_exe
      # Ut_parser_tstParallel_File_Token_Stream_exe
      # Ut_parser_tstParse_Table_exe
      # Ut_parser_tstString_Token_Stream_exe
      # Ut_parser_tstToken_Equivalence_exe
      # Ut_parser_tstToken_exe
      # Ut_parser_tstUnit_exe
      # Ut_parser_tstutilities_exe
      # Ut_roots_tstbroydn_exe
      # Ut_roots_tstcubic1_exe
      # Ut_roots_tstfdjac_exe
      # Ut_roots_tstlnsrch_exe
      # Ut_roots_tstquadratic_exe
      # Ut_roots_tstzbrac_exe
      # Ut_roots_tstzbrent_exe
      # Ut_timestep_tstTimeStep_exe
      # Ut_viz_install_inputs_1
      # Ut_viz_install_inputs_2
      # Ut_viz_install_inputs_3
      # Ut_viz_tstEnsight_Stream_exe
      # Ut_viz_tstEnsight_Translator_exe
      # Ut_cdi_ipcress_ReadOdfIpcressFile_exe
      # Ut_cdi_ipcress_tIpcressFile_exe
      # Ut_cdi_ipcress_tIpcressOdfmgOpacity_exe
      # Ut_cdi_ipcress_tIpcressOpacity_exe
      # Ut_cdi_ipcress_tIpcressWithCDI_exe
      # Ut_cdi_ipcress_tIpcress_Interpreter_exe
      # Ut_special_functions_tstExpInt_exe
      # Ut_special_functions_tstF12_exe
      # Ut_special_functions_tstF12inv_exe
      # Ut_special_functions_tstF1_exe
      # Ut_special_functions_tstF2inv_exe
      # Ut_special_functions_tstF32_exe
      # Ut_special_functions_tstF3_exe
      # Ut_special_functions_tstF4_exe
      # Ut_special_functions_tstFM12_exe
      # Ut_special_functions_tstF_eta_inv_exe
      # Ut_special_functions_tstPower_exe
      # Ut_special_functions_tstSpecialFunctions_exe
      # Ut_special_functions_tstYlm_exe
      # Ut_cdi_analytic_tstAnalytic_EoS_exe
      # Ut_cdi_analytic_tstAnalytic_Gray_Opacity_exe
      # Ut_cdi_analytic_tstAnalytic_Odfmg_Opacity_exe
      # Ut_cdi_analytic_tstPseudo_Line_Analytic_MultigroupOpacity_exe
      # Ut_cdi_analytic_tstnGray_Analytic_MultigroupOpacity_exe
      # Ut_RTT_Format_Reader_TestRTTFormatReader_exe
      # Ut_RTT_Format_Reader_TestRTTMeshReader_exe
      # Ut_RTT_Format_Reader_install_inputs_1
      # Ut_quadrature_tAngle_Operator_exe
      # Ut_quadrature_tAxialQuadrature_exe
      # Ut_quadrature_tQuadCreator_exe
      # Ut_quadrature_tQuadServices_exe
      # Ut_quadrature_tQuadrature_exe
      # Ut_quadrature_tstOrdinate_exe
      # Ut_quadrature_tstgaulag_exe
      # Ut_quadrature_tstgauleg_exe
      Exe_Ipcress_Interpreter
      )
   echo_targets( ${dbs_targets} )
endif()

##---------------------------------------------------------------------------##
## End of src/CMakeLists.txt
##---------------------------------------------------------------------------##
