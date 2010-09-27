# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( fpe_trap_test CXX )

# ---------------------------------------------------------------------------- #
# Prereq
# ---------------------------------------------------------------------------- #
# This unit test requires Python
find_package(PythonInterp)
#  PYTHONINTERP_FOUND - Was the Python executable found
#  PYTHON_EXECUTABLE  - path to the Python interpreter
if( PYTHONINTERP_FOUND )

    # ---------------------------------------------------------------------------- #
    # Source files
    # ---------------------------------------------------------------------------- #

    # set( test_lib_sources 
       # ${PROJECT_SOURCE_DIR}/fpe_trap_test.cc 
       # ${PROJECT_SOURCE_DIR}/fpe_trap_test.hh
    # )
    file( GLOB test_sources *.cc )
    # list( REMOVE_ITEM test_sources ${test_lib_sources} )

    # ---------------------------------------------------------------------------- #
    # Directories to search for include directives
    # ---------------------------------------------------------------------------- #

    include_directories( 
       ${PROJECT_SOURCE_DIR}      # headers for tests
       ${PROJECT_SOURCE_DIR}/..   # headers for package
    )

    # ---------------------------------------------------------------------------- #
    # Build library for test directory
    # ---------------------------------------------------------------------------- #

    # add_library( Lib_fpe_trap_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
    # if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
       # set_target_properties( 
          # Lib_fpe_trap_test
            # PROPERTIES 
            # # Provide compile define macro to enable declspec(dllexport) linkage.
            # COMPILE_DEFINITIONS BUILDING_DLL 
            # # Use custom library naming
            # OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}fpe_trap_test${CMAKE_SHARED_LIBRARY_SUFFIX}
       # )
    # else()
       # set_target_properties( 
          # Lib_fpe_trap_test 
          # PROPERTIES 
            # # Use custom library naming
        # OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}fpe_trap_test${CMAKE_SHARED_LIBRARY_SUFFIX}
       # )
    # endif()

    # For Win32 with shared libraries, the package dll must be located in
    # the test directory.

    # get_target_property( Lib_fpe_trap_loc Lib_fpe_trap LOCATION )
    # if( WIN32 )
       # add_custom_command( TARGET Lib_fpe_trap_test
          # POST_BUILD
          # COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_fpe_trap_loc} 
          # ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
       # )
    # endif()
    # add_dependencies( Lib_fpe_trap_test Lib_fpe_trap )

    # ---------------------------------------------------------------------------- #
    # Build Unit tests
    # ---------------------------------------------------------------------------- #

    add_executable( Ut_fpe_trap_do_exception_exe ${test_sources} )
    set_target_properties( 
       Ut_fpe_trap_do_exception_exe 
       PROPERTIES 
           OUTPUT_NAME do_exception 
           VS_KEYWORD do_exception
           PROJECT_LABEL Ut_fpe_trap
    )
    target_link_libraries( Ut_fpe_trap_do_exception_exe 
      Lib_fpe_trap 
      Lib_dsxx 
    )
    add_test( 
        NAME fpe_trap_python_test 
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/python/test_fpe_trap.py $<TARGET_FILE:Ut_fpe_trap_do_exception_exe>
    )
    set_tests_properties( fpe_trap_python_test 
       PROPERTIES	
         PASS_REGULAR_EXPRESSION ".*PASSED" 
         FAIL_REGULAR_EXPRESSION ".*FAILED" 
    )


endif( PYTHONINTERP_FOUND )
    
# source_group( fpe_trap\\test FILES ${test_lib_sources} ${test_sources} )
