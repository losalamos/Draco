#-----------------------------*-cmake-*----------------------------------------#
# Generic CMake Variable Logging
#
# See: http://www.cmake.org/Wiki/CMake_Useful_Variables
# -----------------------------------------------------------------------------#

MESSAGE( "\n------------------------------------------------------------" )
MESSAGE( "Begin varlogger" )
MESSAGE( "------------------------------------------------------------" )

#########################################
# PATHS
#########################################

MESSAGE( "\nPATHS\n------------------------------------------------------------")

# if you are building in-source, this is the same as CMAKE_SOURCE_DIR,
# otherwise this is the top level directory of your build tree  
MESSAGE( STATUS "CMAKE_BINARY_DIR:         " ${CMAKE_BINARY_DIR} )

# if you are building in-source, this is the same as
# CMAKE_CURRENT_SOURCE_DIR, otherwise this is the directory where the
# compiled or generated files from the current CMakeLists.txt will go
# to  
MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )

# this is the directory, from which cmake was started, i.e. the top
# level source directory  
MESSAGE( STATUS "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR} )

# this is the directory where the currently processed CMakeLists.txt
# is located in  
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

# The name of the project
message( STATUS "PROJECT_NAME            : " ${PROJECT_NAME} )

# contains the full path to the top level directory of your build tree
MESSAGE( STATUS "PROJECT_BINARY_DIR      : " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the
# PROJECT() command.
MESSAGE( STATUS "PROJECT_SOURCE_DIR:       " ${PROJECT_SOURCE_DIR} )

# set this variable to specify a common place where CMake should put
# all executable files (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH:   " ${EXECUTABLE_OUTPUT_PATH} )

# set this variable to specify a common place where CMake should put
# all libraries (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE( STATUS "LIBRARY_OUTPUT_PATH:      " ${LIBRARY_OUTPUT_PATH} )

# tell CMake to search first in directories listed in CMAKE_MODULE_PATH
# when you use FIND_PACKAGE() or INCLUDE()
MESSAGE( STATUS "CMAKE_MODULE_PATH:        " ${CMAKE_MODULE_PATH} )

MESSAGE( STATUS "CMAKE_FRAMEWORK_PATH       : " ${CMAKE_FRAMEWORK_PATH} )
MESSAGE( STATUS "CMAKE_APPBUNDLE_PATH       : " ${CMAKE_APPBUNDLE_PATH} )
MESSAGE( STATUS "CMAKE_INCLUDE_PATH         : " ${CMAKE_INCLUDE_PATH} )
MESSAGE( STATUS "CMAKE_LIBRARY_PATH         : " ${CMAKE_LIBRARY_PATH} )

MESSAGE( STATUS "CMAKE_SYSTEM_FRAMEWORK_PATH: " ${CMAKE_SYSTEM_FRAMEWORK_PATH} )
MESSAGE( STATUS "CMAKE_SYSTEM_APPBUNDLE_PATH: " ${CMAKE_SYSTEM_APPBUNDLE_PATH} )
MESSAGE( STATUS "CMAKE_SYSTEM_INCLUDE_PATH  : " ${CMAKE_SYSTEM_INCLUDE_PATH} )
MESSAGE( STATUS "CMAKE_SYSTEM_LIBRARY_PATH  : " ${CMAKE_SYSTEM_LIBRARY_PATH} )

MESSAGE( STATUS "ENV(PATH)                  :     " $ENV{PATH} )
MESSAGE( STATUS "ENV(LIB)                   :     " $ENV{LIB} )
MESSAGE( STATUS "ENV(INCLUDE)               :     " $ENV{INCLUDE} )


#########################################
# CMAKE Informatoin
#########################################

MESSAGE( "\nCMake Information:\n------------------------------------------------------------")

# this is the complete path of the cmake which runs currently
# (e.g. /usr/local/bin/cmake)  
MESSAGE( STATUS "CMAKE_COMMAND: " ${CMAKE_COMMAND} )

# this is the CMake installation directory 
MESSAGE( STATUS "CMAKE_ROOT: " ${CMAKE_ROOT} )

# this is the filename including the complete path of the file where
# this variable is used.  
MESSAGE( STATUS "CMAKE_CURRENT_LIST_FILE: " ${CMAKE_CURRENT_LIST_FILE} )

# this is linenumber where the variable is used
MESSAGE( STATUS "CMAKE_CURRENT_LIST_LINE: " ${CMAKE_CURRENT_LIST_LINE} )

# this is used when searching for include files e.g. using the
# FIND_PATH() command. 
MESSAGE( STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH} )

# this is used when searching for libraries e.g. using the
# FIND_LIBRARY() command. 
MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

# the directory within the current binary directory that contains all
# the CMake generated files. Typically evaluates to
# "/CMakeFiles". Note the leading slash for the directory. Typically
# used with the current binary directory,
# i.e. ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY} 
MESSAGE( STATUS "CMAKE_FILES_DIRECTORY: "${CMAKE_FILES_DIRECTORY} )

message( STATUS "CMAKE_MAKE_PROGRAM:    ${CMAKE_MAKE_PROGRAM}" )


########################################
# System Information
########################################

MESSAGE( "\nSystem Information:\n------------------------------------------------------------")

# the complete system name, e.g. "Linux-2.6.9-42.0.10.ELsmp"
# or "Windows 5.1"  
MESSAGE( STATUS "CMAKE_SYSTEM: " ${CMAKE_SYSTEM} )

# the short system name, e.g. "Linux", "FreeBSD" or "Windows"
MESSAGE( STATUS "CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME} )

# only the version part of CMAKE_SYSTEM, e.g. "2.6.9-42.0.10.ELsmp"
MESSAGE( STATUS "CMAKE_SYSTEM_VERSION: " ${CMAKE_SYSTEM_VERSION} )

# the processor name, e.g. "x86_64," "i686"
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR: " ${CMAKE_SYSTEM_PROCESSOR} )

# is TRUE on all UNIX-like OS's, including Apple OS X and CygWin
MESSAGE( STATUS "UNIX: " ${UNIX} )

# is TRUE on Windows, including CygWin 
MESSAGE( STATUS "WIN32: " ${WIN32} )

# is TRUE on Apple OS X
MESSAGE( STATUS "APPLE: " ${APPLE} )

# is TRUE when using the MinGW compiler in Windows
MESSAGE( STATUS "MINGW: " ${MINGW} )

# is TRUE on Windows when using the CygWin version of cmake
MESSAGE( STATUS "CYGWIN: " ${CYGWIN} )

# is TRUE on Windows when using a Borland compiler 
MESSAGE( STATUS "BORLAND: " ${BORLAND} )

########################################
# Compiler Information
########################################

MESSAGE( "\nCompiler Information:\n------------------------------------------------------------")

# Microsoft compiler 
MESSAGE( STATUS "MSVC: " ${MSVC} )
MESSAGE( STATUS "MSVC_IDE: " ${MSVC_IDE} )
MESSAGE( STATUS "MSVC60: " ${MSVC60} )
MESSAGE( STATUS "MSVC70: " ${MSVC70} )
MESSAGE( STATUS "MSVC71: " ${MSVC71} )
MESSAGE( STATUS "MSVC80: " ${MSVC80} )
MESSAGE( STATUS "CMAKE_COMPILER_2005: " ${CMAKE_COMPILER_2005} )

# set this to true if you don't want to rebuild the object files if
# the rules have changed, but not the actual source files or headers
# (e.g. if you changed the some compiler switches)
MESSAGE( STATUS "CMAKE_SKIP_RULE_DEPENDENCY: " 
  ${CMAKE_SKIP_RULE_DEPENDENCY} )

# since CMake 2.1 the install rule depends on all, i.e. everything
# will be built before installing. If you don't like this, set this
# one to true.
MESSAGE( STATUS "CMAKE_SKIP_INSTALL_ALL_DEPENDENCY: " 
  ${CMAKE_SKIP_INSTALL_ALL_DEPENDENCY} )

# If set, runtime paths are not added when using shared
# libraries. Default it is set to OFF
MESSAGE( STATUS "CMAKE_SKIP_RPATH: " ${CMAKE_SKIP_RPATH} )

# set this to true if you are using makefiles and want to see the full
# compile and link commands instead of only the shortened ones
MESSAGE( STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )

# this will cause CMake to not put in the rules that re-run
# CMake. This might be useful if you want to use the generated build
# files on another machine.
MESSAGE( STATUS "CMAKE_SUPPRESS_REGENERATION: " 
  ${CMAKE_SUPPRESS_REGENERATION} )

# A simple way to get switches to the compiler is to use
# ADD_DEFINITIONS(). But there are also two variables exactly for this
# purpose:

# the compiler flags for compiling C sources 
MESSAGE( STATUS "CMAKE_C_FLAGS: " ${CMAKE_C_FLAGS} )

# the compiler flags for compiling C++ sources 
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )


# Choose the type of build.  Example: SET(CMAKE_BUILD_TYPE Debug)
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

# if this is set to ON, then all libraries are built as shared
# libraries by default.
MESSAGE( STATUS "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS} )

# the compiler used for C files 
MESSAGE( STATUS "CMAKE_C_COMPILER: " ${CMAKE_C_COMPILER} )

# the compiler used for C++ files 
MESSAGE( STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )

# if the compiler is a variant of gcc, this should be set to 1 
MESSAGE( STATUS "CMAKE_COMPILER_IS_GNUCC: " ${CMAKE_COMPILER_IS_GNUCC} )

# if the compiler is a variant of g++, this should be set to 1 
MESSAGE( STATUS "CMAKE_COMPILER_IS_GNUCXX : " ${CMAKE_COMPILER_IS_GNUCXX} )

# the tools for creating libraries 
MESSAGE( STATUS "CMAKE_AR: " ${CMAKE_AR} )
MESSAGE( STATUS "CMAKE_RANLIB: " ${CMAKE_RANLIB} )

# Library names
MESSAGE( STATUS "CMAKE_STATIC_LIBRARY_PREFIX: " ${CMAKE_STATIC_LIBRARY_PREFIX} )
MESSAGE( STATUS "CMAKE_STATIC_LIBRARY_SUFFIX: " ${CMAKE_STATIC_LIBRARY_SUFFIX} )
MESSAGE( STATUS "CMAKE_SHARED_LIBRARY_PREFIX: " ${CMAKE_SHARED_LIBRARY_PREFIX} )
MESSAGE( STATUS "CMAKE_SHARED_LIBRARY_SUFFIX: " ${CMAKE_SHARED_LIBRARY_SUFFIX} )
MESSAGE( STATUS "CMAKE_SHARED_MODULE_PREFIX : " ${CMAKE_SHARED_MODULE_PREFIX} )
MESSAGE( STATUS "CMAKE_SHARED_MODULE_SUFFIX : " ${CMAKE_SHARED_MODULE_SUFFIX} )
message( STATUS "CMAKE_EXECUTABLE_SUFFIX    : " ${CMAKE_EXECUTABLE_SUFFIX} )

#
#MESSAGE( STATUS ": " ${} )

# CFLAGS
message("CMAKE_DL_LIBS = ${CMAKE_DL_LIBS}")
message("CMAKE_SHARED_LIBRARY_C_FLAGS = ${CMAKE_SHARED_LIBRARY_C_FLAGS}")
message("CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS = ${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
message("CMAKE_SHARED_LIBRARY_LINK_C_FLAGS  = ${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS}")
message("CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG = ${CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG}")
message("CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP = ${CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP}")
message("CMAKE_SHARED_LIBRARY_RPATH_LINK_C_FLAG = ${CMAKE_SHARED_LIBRARY_RPATH_LINK_C_FLAG}")
message("CMAKE_SHARED_LIBRARY_SONAME_C_FLAG = ${CMAKE_SHARED_LIBRARY_SONAME_C_FLAG}")
message("CMAKE_SHARED_LIBRARY_SONAME_CXX_FLAG = ${CMAKE_SHARED_LIBRARY_SONAME_CXX_FLAG}")
message("CMAKE_SHARED_LIBRARY_SONAME_Fortran_FLAG  = ${CMAKE_SHARED_LIBRARY_SONAME_Fortran_FLAG}")
message("CMAKE_EXE_EXPORTS_C_FLAG  = ${CMAKE_EXE_EXPORTS_C_FLAG}")
message("CMAKE_EXE_EXPORTS_CXX_FLAG  = ${CMAKE_EXE_EXPORTS_CXX_FLAG}")

MESSAGE( "\nFortran Compiler Information:\n------------------------------------------------------------")
# message("CMAKE_Fortran_FLAGS_INIT                = ${CMAKE_Fortran_FLAGS_INIT}" )
# message("CMAKE_Fortran_FLAGS_DEBUG_INIT          = ${CMAKE_Fortran_FLAGS_DEBUG_INIT}" )
# message("CMAKE_Fortran_FLAGS_RELEASE_INIT        = ${CMAKE_Fortran_FLAGS_RELEASE_INIT}" )
# message("CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT = ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT}" )
message("CMAKE_Fortran_FLAGS                = ${CMAKE_Fortran_FLAGS}" )
message("CMAKE_Fortran_FLAGS_DEBUG          = ${CMAKE_Fortran_FLAGS_DEBUG}" )
message("CMAKE_Fortran_FLAGS_RELEASE        = ${CMAKE_Fortran_FLAGS_RELEASE}" )
message("CMAKE_Fortran_FLAGS_RELWITHDEBINFO = ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" )
message("CMAKE_Fortran_FLAGS_MINSIZEREL     = ${CMAKE_Fortran_FLAGS_MINSIZEREL}" )

message("CMAKE_Fortran_STANDARD_LIBRARIES          = ${MAKE_Fortran_STANDARD_LIBRARIES}" )
message("CMAKE_SHARED_LIBRARY_Fortran_FLAGS        = ${CMAKE_SHARED_LIBRARY_Fortran_FLAGS}")
message("CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS = ${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
message("CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS   = ${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")

message("CMAKE_Fortran_CREATE_SHARED_LIBRARY = ${CMAKE_Fortran_CREATE_SHARED_LIBRARY}")
message("CMAKE_Fortran_CREATE_STATIC_LIBRARY = ${CMAKE_Fortran_CREATE_STATIC_LIBRARY}")
message("CMAKE_Fortran_COMPILE_OBJECT        = ${CMAKE_Fortran_COMPILE_OBJECT}")
message("CMAKE_COMPILE_RESOURCE              = ${CMAKE_COMPILE_RESOURCE}")
message("CMAKE_Fortran_LINK_EXECUTABLE       = ${CMAKE_Fortran_LINK_EXECUTABLE}")


# EXE_LINKER
message("CMAKE_LINKER           = ${CMAKE_LINKER}" )
message("CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
message("CMAKE_EXE_LINKER_FLAGS_DEBUG          = ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message("CMAKE_EXE_LINKER_FLAGS_RELEASE        = ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
message("CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO = ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}")
message("CMAKE_EXE_LINKER_FLAGS_MINSIZEREL     = ${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}")

message("CMAKE_MODULE_LINKER_FLAGS         = ${CMAKE_MODULE_LINKER_FLAGS}")
message("CMAKE_MODULE_LINKER_FLAGS_DEBUG   = ${CMAKE_MODULE_LINKER_FLAGS_DEBUG}")
message("CMAKE_MODULE_LINKER_FLAGS_RELEASE = ${CMAKE_MODULE_LINKER_FLAGS_RELEASE}")
message("CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO = ${CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO}")
message("CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL     = ${CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL}")

message("CMAKE_SHARED_LINKER_FLAGS         = ${CMAKE_SHARED_LINKER_FLAGS}")
message("CMAKE_SHARED_LINKER_FLAGS_DEBUG   = ${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")
message("CMAKE_SHARED_LINKER_FLAGS_RELEASE = ${CMAKE_SHARED_LINKER_FLAGS_RELEASE}")
message("CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO = ${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO}")
message("CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL     = ${CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL}")


########################################
# Data Types and Sizes
########################################
# CMakeBackwardCompatibilityC.cmake

MESSAGE( "\nData Types and Sizes:\n------------------------------------------------------------")

# Size of default types
#MESSAGE( STATUS "CMAKE_SIZEOF_INT    : " ${CMAKE_SIZEOF_INT} )
#MESSAGE( STATUS "CMAKE_SIZEOF_LONG   : " ${CMAKE_SIZEOF_LONG} )
MESSAGE( STATUS "CMAKE_SIZEOF_VOID_P : " ${CMAKE_SIZEOF_VOID_P} )
#MESSAGE( STATUS "CMAKE_SIZEOF_CHAR   : " ${CMAKE_SIZEOF_CHAR} )
#MESSAGE( STATUS "CMAKE_SIZEOF_SHORT  : " ${CMAKE_SIZEOF_SHORT} )
#MESSAGE( STATUS "CMAKE_SIZEOF_FLOAT  : " ${CMAKE_SIZEOF_FLOAT} )
#MESSAGE( STATUS "CMAKE_SIZEOF_DOUBLE : " ${CMAKE_SIZEOF_DOUBLE} )

# ------------------------- End of Generic CMake Variable Logging ------------------
