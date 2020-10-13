#-----------------------------*-cmake-*----------------------------------------#
# file   compton2/test/CMakeLists.txt
# author Andrew Till
# date   11 May 2020
# brief  Generate build project files for compton2/test.
# note   Copyright (C) 2017-2020, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

project( compton2_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
  ${PROJECT_SOURCE_DIR}/tCompton2.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_parallel_tests(
    SOURCES "${test_sources}"
    PE_LIST "1"
    DEPS    "Lib_compton2;Lib_c4"
    LABEL   "nomemcheck"
    MPI_PLUS_OMP )

#add_app_unit_test(
#    DRIVER    ${CMAKE_CURRENT_SOURCE_DIR}/tIpcress_Interpreter.py
#    APP       $<TARGET_FILE_DIR:Exe_Ipcress_Interpreter>/$<TARGET_FILE_NAME:Exe_Ipcress_Interpreter>
#    TEST_ARGS "${PROJECT_SOURCE_DIR}/odfregression10.ipcress"
#    STDINFILE "${PROJECT_SOURCE_DIR}/IpcressInterpreter.odf.stdin"
#    GOLDFILE  "${PROJECT_SOURCE_DIR}/Ipcress_Interpreter-odfregression10.ipcress.gold"
#    )
#add_test(
#    NAME Ut_compton2_cskrw
#    COMMAND $<TARGET_FILE:CSKrw> dummy_data
#    )

# ---------------------------------------------------------------------------- #
# end of compton2/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
