# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( lapack_wrap_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/lapack_wrap_test.cc 
   ${PROJECT_SOURCE_DIR}/lapack_wrap_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

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

add_component_library( Lib_lapack_wrap_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies(      Lib_lapack_wrap_test Lib_lapack_wrap )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_dsxx;${LAPACK_LIBRARIES};${FORTRAN_LIBS}" )

# foreach( file ${test_sources} )
#   get_filename_component( testname ${file} NAME_WE )
#   add_executable( Ut_lapack_wrap_${testname}_exe ${file})
#   set_target_properties( 
#      Ut_lapack_wrap_${testname}_exe 
#      PROPERTIES 
#        OUTPUT_NAME ${testname} 
#        VS_KEYWORD ${testname}
#        PROJECT_LABEL Ut_lapack_wrap
#        #LINK_FLAGS /nodefaultlib:msvcrt.lib
#   )
#   target_link_libraries( Ut_lapack_wrap_${testname}_exe Lib_lapack_wrap Lib_lapack_wrap_test Lib_dsxx ${LAPACK_LIBRARIES} ${FORTRAN_LIBS} )
#   add_test( lapack_wrap_${testname} ${testname} )
#   set_tests_properties( lapack_wrap_${testname} PROPERTIES	
#      PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED" 
#      FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
#   )
#   if( "${test_loc}x" STREQUAL "x" )
#     get_target_property( test_loc Ut_lapack_wrap_${testname}_exe LOCATION )
#     get_filename_component( test_loc ${test_loc} PATH )
#   endif()
# endforeach()


