#-----------------------------*-cmake-*----------------------------------------#
# file   viz/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for viz/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( viz_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources
  ${PROJECT_SOURCE_DIR}/tstEnsight_Stream.cc
  ${PROJECT_SOURCE_DIR}/tstViz_Traits.cc
  )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

if( NOT APPLE )
  # numdiff/ndselect crashes when run on OS/X so don't try to run the tests.

  add_scalar_tests(
    SOURCES "${test_sources}"
    DEPS    "Lib_viz" )

  # for Ut_viz_tstEnsight_Translator_exe, check the generated viz files using
  # numdiff from a cmake script instead of trying to use the system() command.
  # Using system fails for BG/Q systems.

  add_component_executable(
    TARGET      Ut_viz_tstEnsight_Translator_exe
    TARGET_DEPS Lib_viz
    EXE_NAME    tstEnsight_Translator
    SOURCES     ${PROJECT_SOURCE_DIR}/tstEnsight_Translator.cc
    PREFIX      Draco )

  add_component_executable(
    TARGET      Ut_viz_tstEnsight_Translator_Unstructured_exe
    TARGET_DEPS Lib_viz
    EXE_NAME    tstEnsight_Translator_Unstructured
    SOURCES     ${PROJECT_SOURCE_DIR}/tstEnsight_Translator_Unstructured.cc
    PREFIX      Draco )

  # set( VERBOSE_DEBUG ON )
  include( ApplicationUnitTest )
  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tstEnsight_Translator.py
    APP    $<TARGET_FILE_DIR:Ut_viz_tstEnsight_Translator_exe>/$<TARGET_FILE_NAME:Ut_viz_tstEnsight_Translator_exe>
    BUILDENV "CMAKE_HOST_SYSTEM_PROCESSOR=${CMAKE_HOST_SYSTEM_PROCESSOR}"
    )

  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tstEnsight_Translator_Unstructured.py
    APP    $<TARGET_FILE_DIR:Ut_viz_tstEnsight_Translator_Unstructured_exe>/$<TARGET_FILE_NAME:Ut_viz_tstEnsight_Translator_Unstructured_exe>
    BUILDENV "CMAKE_HOST_SYSTEM_PROCESSOR=${CMAKE_HOST_SYSTEM_PROCESSOR}"
    )

endif()

##---------------------------------------------------------------------------##
## Extra 'make clean' files
##---------------------------------------------------------------------------##
set( extra_clean_files
   Densities.diff
   Pressure.diff
   Temperatures.diff
   Velocity.diff
   ensight_stream.out
   geo.diff
   part_testproblem_binary_ensight
   part_testproblem_ensight
   testproblem_binary_ensight
   testproblem_ensight
)
set_directory_properties(
   PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${extra_clean_files}" )

#------------------------------------------------------------------------------#
# End viz/test/CMakeLists.txt
#------------------------------------------------------------------------------#
