##---------------------------------------------------------------------------##
# file   : platform_checks.cmake
# brief  : Platform Checks for Draco Build System
# note   : Copyright (C) 2016 Los Alamos National Security, LLC.
#          All rights reserved
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Determine System Type and System Names
##
## Used by ds++ and c4.
##---------------------------------------------------------------------------##
macro( set_draco_uname )
    # Store platform information in config.h
    if( UNIX )
      set( draco_isLinux 1 )
      set( DRACO_UNAME "Linux" )
    elseif( WIN32 )
      set( draco_isWin 1 )
      set( DRACO_UNAME "Windows" )
    elseif( OSF1 )
      set( draco_isOSF1 1 )
      set( DRACO_UNAME "OSF1" )
    elseif( APPLE )
      set( draco_isDarwin 1)
      set( DRACO_UNAME "Darwin" )
    else()
      set( draco_isAIX 1 )
      set( DRACO_UNAME "AIX" )
    endif()
    # Special setup for catamount
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Catamount" )
       set( draco_isLinux_with_aprun 1 )
       set( draco_isCatamount 1 )
       set( DRACO_UNAME "Catamount" )
    endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine if gethostname() is available.
## Determine the value of HOST_NAME_MAX.
## 
## Used by ds++/SystemCall.cc and ds++/path.cc
##---------------------------------------------------------------------------##
macro( query_have_gethostname )
    # Platform checks for gethostname()
    include( CheckIncludeFiles )
    check_include_files( unistd.h    HAVE_UNISTD_H )
    check_include_files( limits.h    HAVE_LIMITS_H )
    check_include_files( winsock2.h  HAVE_WINSOCK2_H )
    check_include_files( direct.h    HAVE_DIRECT_H )
    check_include_files( sys/param.h HAVE_SYS_PARAM_H )
    # Used to demangle symbols for stack trace
    # check_include_files( cxxabi.h    HAVE_CXXABI_H )

    # -------------- Checks for hostname and len(hostname) ---------------- #
    # gethostname()
    include( CheckFunctionExists )
    check_function_exists( gethostname HAVE_GETHOSTNAME )

    # HOST_NAME_MAX
    include( CheckSymbolExists )
    unset( hlist )
    if( HAVE_UNISTD_H )
       list( APPEND hlist unistd.h )
    endif()
    if( HAVE_WINSOCK2_H )
       list( APPEND hlist winsock2.h )
    endif()
    if( HAVE_LIMITS_H )
       list( APPEND hlist limits.h )
    endif()
    check_symbol_exists( HOST_NAME_MAX "${hlist}" HAVE_HOST_NAME_MAX )
    if( NOT HAVE_HOST_NAME_MAX )
       unset( HAVE_GETHOSTNAME )
    endif()

    check_symbol_exists( _POSIX_HOST_NAME_MAX "posix1_lim.h" HAVE_POSIX_HOST_NAME_MAX )

    # HOST_NAME_MAX
    check_symbol_exists( MAXHOSTNAMELEN "sys/param.h" HAVE_MAXHOSTNAMELEN )
    if( NOT HAVE_MAXHOSTNAMELEN )
       unset( HAVE_MAXHOSTNAMELEN )
    endif()
    
endmacro()

##---------------------------------------------------------------------------##
## Determine if gethostname() is available.
## Determine the value of HOST_NAME_MAX.
## 
## Used by ds++/SystemCall.cc and ds++/path.cc
##---------------------------------------------------------------------------##
macro( query_have_maxpathlen )
    # MAXPATHLEN
    unset( hlist )
    if( HAVE_UNISTD_H )
       list( APPEND hlist unistd.h )
    endif()
    if( HAVE_LIMITS_H )
       list( APPEND hlist limits.h )
    endif()
    if( HAVE_SYS_PARAM_H )
       list( APPEND hlist sys/param.h )
    endif()
    check_symbol_exists( MAXPATHLEN "${hlist}" HAVE_MAXPATHLEN )
    if( NOT HAVE_MAXPATHLEN )
        unset( HAVE_MAXPATHLEN )
    endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine if some system headers exist
##---------------------------------------------------------------------------##
macro( query_have_sys_headers )
   
   include( CheckIncludeFiles )
   check_include_files( sys/types.h HAVE_SYS_TYPES_H )
   check_include_files( unistd.h    HAVE_UNISTD_H    )

endmacro()

##---------------------------------------------------------------------------##
## Check 8-byte int type
##
## For some systems, provide special compile flags to support 8-byte integers
##---------------------------------------------------------------------------##

macro(check_eight_byte_int_type)
   if( "${SIZEOF_INT}notset" STREQUAL "notset" )
      determine_word_types()
   endif()

   if( "${SIZEOF_INT}" STREQUAL "8" )
      message( "Checking for 8-byte integer type... int - no mods needed." )
   elseif( "${SIZEOF_LONG}" STREQUAL "8" )
      message( "Checking for 8-byte integer type... long - no mods needed." )
   else()
      message( FATAL_ERROR "need to patch up this part of the build system." )
   endif()
endmacro()

##---------------------------------------------------------------------------##
## Detect support for the C99 restrict keyword
## Borrowed from http://cmake.3232098.n2.nabble.com/AC-C-RESTRICT-td7582761.html
##
## A restrict-qualified pointer (or reference) is basically a promise to the 
## compiler that for the scope of the pointer, the target of the pointer will 
## only be accessed through that pointer (and pointers copied from it).
##
## http://www.research.scea.com/research/pdfs/GDC2003_Memory_Optimization_18Mar03.pdf
##---------------------------------------------------------------------------##
macro( query_have_restrict_keyword )
   
   message(STATUS "Looking for the C99 restrict keyword")
   include( CheckCSourceCompiles )
   foreach( ac_kw __restrict __restrict__ _Restrict restrict )
      check_c_source_compiles("
         typedef int * int_ptr;
         int foo ( int_ptr ${ac_kw} ip ) { return ip[0]; }
         int main (void) {
            int s[1];
            int * ${ac_kw} t = s;
            t[0] = 0;
            return foo(t); }
         "
         HAVE_RESTRICT) 

      if( HAVE_RESTRICT )
         set( RESTRICT_KEYWORD ${ac_kw} )
         message(STATUS "Looking for the C99 restrict keyword - found ${RESTRICT_KEYWORD}")
         break()
      endif()
   endforeach()
   if( NOT HAVE_RESTRICT )
      message(STATUS "Looking for the C99 restrict keyword - not found")
   endif()

endmacro()

##---------------------------------------------------------------------------##
## Detect C++11 features
##
## This macro requires CMake 3.1+
##
## 1. This macro detects available C++11 features and sets CPP macros in
##    the build system that have the form HAS_CXX11_<FEATURE>.  These
##    values are saved in ds++/config.h.
## 2. This macro also checks to ensure that the current C++ compiler
##    supports the C++11 features already in use by draco.
##
## http://stackoverflow.com/questions/23042722/how-to-detect-which-c11-features-are-used-in-my-source-code
## http://www.cmake.org/cmake/help/v3.1/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html
##
## CMake will automatically add the '-std=c++11' compiler flag if it
## sees a command of the following form:
##
## target_compile_features( Lib_dsxx PRIVATE cxx_auto_type )
##
## Draco adds this flag automatically in config/unix-g++.cmake so the
## above probably isn't needed anywhere.
##
##---------------------------------------------------------------------------##
macro( query_cxx11_features )

  message( STATUS "Looking for required C++11 features..." )
  get_property(cxx_features GLOBAL PROPERTY CMAKE_CXX_KNOWN_FEATURES)
  # compatibility with the old C++11 feature detection system
  set( CXX11_FEATURE_LIST "${cxx_features}" CACHE STRING "List of known C++11 features (ds++/config.h)." FORCE )
  set( cxx11_required_features
    cxx_auto_type
    cxx_decltype_auto
#    cxx_nullptr
#    cxx_lambdas
    cxx_rvalue_references
    cxx_long_long_type
    cxx_static_assert
    cxx_decltype
#    cxx_variadic_templates
#    cxx_sizeof_member
#    cxx_generalized_initializers
    )
# cxx_aggregate_default_initializers
# cxx_alias_templates
# cxx_alignas
# cxx_alignof
# cxx_attributes
# cxx_attribute_deprecated
# cxx_auto_type
# cxx_binary_literals
# cxx_constexpr
# cxx_contextual_conversions
# cxx_decltype
# cxx_decltype_auto
# cxx_decltype_incomplete_return_types
# cxx_default_function_template_args
# cxx_defaulted_functions
# cxx_defaulted_move_initializers
# cxx_delegating_constructors
# cxx_deleted_functions
# cxx_digit_separators
# cxx_enum_forward_declarations
# cxx_explicit_conversions
# cxx_extended_friend_declarations
# cxx_extern_templates
# cxx_final
# cxx_func_identifier
# cxx_generalized_initializers
# cxx_generic_lambdas
# cxx_inheriting_constructors
# cxx_inline_namespaces
# cxx_lambdas
# cxx_lambda_init_captures
# cxx_local_type_template_args
# cxx_long_long_type
# cxx_noexcept
# cxx_nonstatic_member_init
# cxx_nullptr
# cxx_override
# cxx_range_for
# cxx_raw_string_literals
# cxx_reference_qualified_functions
# cxx_relaxed_constexpr
# cxx_return_type_deduction
# cxx_right_angle_brackets
# cxx_rvalue_references
# cxx_sizeof_member
# cxx_static_assert
# cxx_strong_enums
# cxx_template_template_parameters
# cxx_thread_local
# cxx_trailing_return_types
# cxx_unicode_literals
# cxx_uniform_initialization
# cxx_unrestricted_unions
# cxx_user_literals
# cxx_variable_templates
# cxx_variadic_macros
# cxx_variadic_templates
    
  foreach( cxx11reqfeature ${cxx11_required_features} )
    string( TOUPPER ${cxx11reqfeature} reqfeat )
    string( REPLACE "CXX_" "HAS_CXX11_" reqfeat ${reqfeat} )
    if( NOT "${cxx_features}" MATCHES "${cxx11reqfeature}" )
      message( FATAL_ERROR "Draco requires a C++ compiler that can support the '${cxx11reqfeature}' feature of the C++11 standard.")
    endif()
    # if not available, the variable will not be defined
    set( "${reqfeat}" ON CACHE BOOL "C++11 feature macro value." FORCE )
  endforeach()
  message( STATUS "Looking for required C++11 features...done.  See ds++/config.h for details." )
  
endmacro()

##---------------------------------------------------------------------------##
## Query OpenMP availability
##
## This feature is usually compiler specific and a compile flag must be
## added.  For this to work the <platform>-<compiler>.cmake files (eg.
## unix-g++.cmake) call this macro.
##---------------------------------------------------------------------------##
macro( query_openmp_availability )
  message( STATUS "Looking for OpenMP...")
  find_package(OpenMP QUIET)
  if( OPENMP_FOUND )
    message(STATUS "Looking for OpenMP... ${OpenMP_C_FLAGS}")
    set( OPENMP_FOUND ${OPENMP_FOUND} CACHE BOOL "Is OpenMP availalable?" FORCE )
  else()
    message(STATUS "Looking for OpenMP... not found")
  endif()
endmacro()

##---------------------------------------------------------------------------##
## Sample platform checks
##---------------------------------------------------------------------------##

# # Check for nonblocking collectives
# check_function_exists(MPI_Iallgather HAVE_MPI3_NONBLOCKING_COLLECTIVES)
# check_function_exists(MPIX_Iallgather HAVE_MPIX_NONBLOCKING_COLLECTIVES)

# # Check for MPI_IN_PLACE (essentially MPI2 support)
# include(CheckCXXSourceCompiles)
# set(MPI_IN_PLACE_CODE
#     "#include \"mpi.h\"
#      int main( int argc, char* argv[] )
#      {
#          MPI_Init( &argc, &argv );
#          float a;
#          MPI_Allreduce
#          ( MPI_IN_PLACE, &a, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
#          MPI_Finalize();
#          return 0;
#      }
#     ")
# set(CMAKE_REQUIRED_FLAGS ${CXX_FLAGS})
# set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
# set(CMAKE_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})
# check_cxx_source_compiles("${MPI_IN_PLACE_CODE}" HAVE_MPI_IN_PLACE)

# # Look for restrict support
# set(RESTRICT_CODE
#     "int main(void)
#      {
#          int* RESTRICT a;
#          return 0;
#      }")
# set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict__")
# check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict__)
# if(HAVE___restrict__)
# ...
# endif()



##---------------------------------------------------------------------------##



