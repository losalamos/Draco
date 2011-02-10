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

    file( GLOB test_sources *.cc )

    # ---------------------------------------------------------------------------- #
    # Directories to search for include directives
    # ---------------------------------------------------------------------------- #

    include_directories( 
       ${PROJECT_SOURCE_DIR}      # headers for tests
       ${PROJECT_SOURCE_DIR}/..   # headers for package
    )

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
        COMMAND ${PYTHON_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/python/test_fpe_trap.py 
                $<TARGET_FILE:Ut_fpe_trap_do_exception_exe>
    )
    set_tests_properties( fpe_trap_python_test 
       PROPERTIES	
         PASS_REGULAR_EXPRESSION ".*PASSED" 
         FAIL_REGULAR_EXPRESSION ".*FAILED" 
         LABEL                   "nomemcheck"
    )

endif( PYTHONINTERP_FOUND )
    
# source_group( fpe_trap\\test FILES ${test_lib_sources} ${test_sources} )
