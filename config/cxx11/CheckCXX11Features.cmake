# Checks for C++11 features
#  CXX11_FEATURE_LIST - a list containing all supported features
#  HAS_CXX11_AUTO               - auto keyword
#  HAS_CXX11_NULLPTR            - nullptr
#  HAS_CXX11_LAMBDA             - lambdas
#  HAS_CXX11_STATIC_ASSERT      - static_assert()
#  HAS_CXX11_RVALUE_REFERENCES  - rvalue references
#  HAS_CXX11_DECLTYPE           - decltype keyword
#  HAS_CXX11_CSTDINT_H          - cstdint header
#  HAS_CXX11_LONG_LONG          - long long signed & unsigned types
#  HAS_CXX11_VARIADIC_TEMPLATES - variadic templates
#  HAS_CXX11_CONSTEXPR          - constexpr keyword
#  HAS_CXX11_SIZEOF_MEMBER      - sizeof() non-static members
#  HAS_CXX11_FUNC               - __func__ preprocessor constant
#  HAS_CXX11_INITIALIZER_LISTS  - 
# 
# Original script by Rolf Eike Beer
# Modifications by Andreas Weis
#
cmake_minimum_required(VERSION 2.8.3)

macro(CXX11_CHECK_FEATURE FEATURE_NAME FEATURE_NUMBER RESULT_VAR)
   if (NOT DEFINED ${RESULT_VAR})
      set(_bindir "${CMAKE_CURRENT_BINARY_DIR}/cxx11/cxx11_${FEATURE_NAME}")

      if (${FEATURE_NUMBER})
	 set(_SRCFILE_BASE ${CMAKE_CURRENT_LIST_DIR}/c++11-test-${FEATURE_NAME}-N${FEATURE_NUMBER})
	 set(_LOG_NAME "\"${FEATURE_NAME}\" (N${FEATURE_NUMBER})")
      else (${FEATURE_NUMBER})
	 set(_SRCFILE_BASE ${CMAKE_CURRENT_LIST_DIR}/c++11-test-${FEATURE_NAME})
	 set(_LOG_NAME "\"${FEATURE_NAME}\"")
      endif (${FEATURE_NUMBER})
      message(STATUS "Checking C++11 support for ${_LOG_NAME}")

      set(_SRCFILE "${_SRCFILE_BASE}.cc")
      set(_SRCFILE_FAIL "${_SRCFILE_BASE}_fail.cc")
      set(_SRCFILE_FAIL_COMPILE "${_SRCFILE_BASE}_fail_compile.cc")

      if (CMAKE_CROSSCOMPILING)
	 try_compile(${RESULT_VAR} "${_bindir}" "${_SRCFILE}")
	 if (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
	    try_compile(${RESULT_VAR} "${_bindir}_fail" "${_SRCFILE_FAIL}")
	 endif (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
      else (CMAKE_CROSSCOMPILING)
	 try_run(_RUN_RESULT_VAR _COMPILE_RESULT_VAR
	    "${_bindir}" "${_SRCFILE}")
	 if (_COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR)
	    set(${RESULT_VAR} TRUE)
	 else (_COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR)
	    set(${RESULT_VAR} FALSE)
	 endif (_COMPILE_RESULT_VAR AND NOT _RUN_RESULT_VAR)
	 if (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
	    try_run(_RUN_RESULT_VAR _COMPILE_RESULT_VAR
	       "${_bindir}_fail" "${_SRCFILE_FAIL}")
	    if (_COMPILE_RESULT_VAR AND _RUN_RESULT_VAR)
	       set(${RESULT_VAR} TRUE)
	    else (_COMPILE_RESULT_VAR AND _RUN_RESULT_VAR)
	       set(${RESULT_VAR} FALSE)
	    endif (_COMPILE_RESULT_VAR AND _RUN_RESULT_VAR)
	 endif (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL})
      endif (CMAKE_CROSSCOMPILING)
      if (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE})
	 try_compile(_TMP_RESULT "${_bindir}_fail_compile" "${_SRCFILE_FAIL_COMPILE}")
	 if (_TMP_RESULT)
	    set(${RESULT_VAR} FALSE)
	 else (_TMP_RESULT)
	    set(${RESULT_VAR} TRUE)
	 endif (_TMP_RESULT)
      endif (${RESULT_VAR} AND EXISTS ${_SRCFILE_FAIL_COMPILE})

      if (${RESULT_VAR})
	 message(STATUS "Checking C++11 support for ${_LOG_NAME} -- works")
	 LIST(APPEND CXX11_FEATURE_LIST ${RESULT_VAR})
      else (${RESULT_VAR})
	 message(STATUS "Checking C++11 support for ${_LOG_NAME} -- not supported")
      endif (${RESULT_VAR})
      set(${RESULT_VAR} ${${RESULT_VAR}} CACHE INTERNAL "C++11 support for ${_LOG_NAME}")
   endif (NOT DEFINED ${RESULT_VAR})
endmacro(CXX11_CHECK_FEATURE)

# Numbers are paper number (see
# http://wiki.apache.org/stdcxx/C++0xCompilerSupport )

#cxx11_check_feature("auto"                2546 HAS_CXX11_AUTO)
#cxx11_check_feature("nullptr"             2431 HAS_CXX11_NULLPTR)
#cxx11_check_feature("lambda"              2927 HAS_CXX11_LAMBDA)
#cxx11_check_feature("static_assert"       1720 HAS_CXX11_STATIC_ASSERT)
#cxx11_check_feature("rvalue_references"   2118 HAS_CXX11_RVALUE_REFERENCES)
#cxx11_check_feature("decltype"            2343 HAS_CXX11_DECLTYPE)
cxx11_check_feature("cstdint"             ""   HAS_CXX11_CSTDINT_H)
#cxx11_check_feature("long_long"           1811 HAS_CXX11_LONG_LONG)
cxx11_check_feature("variadic_templates"  2555 HAS_CXX11_VARIADIC_TEMPLATES)
cxx11_check_feature("constexpr"           2235 HAS_CXX11_CONSTEXPR)
cxx11_check_feature("sizeof_member"       2253 HAS_CXX11_SIZEOF_MEMBER)
cxx11_check_feature("initializer_lists"   2672 HAS_CXX11_INITIALIZER_LISTS)
cxx11_check_feature("shared_ptr"          ""   HAS_CXX11_SHARED_PTR)
cxx11_check_feature("array"               ""   HAS_CXX11_ARRAY)
cxx11_check_feature("explicit_conversion" 2437 HAS_CXX11_EXPLICIT_CONVERSION)
cxx11_check_feature("unrestricted_unions" 2544 HAS_CXX11_UNRESTRICTED_UNIONS)
cxx11_check_feature("type_traits"         1836 HAS_CXX11_TYPE_TRAITS)

# Please add a CPP macro def in ds++/config.h.

set(CXX11_FEATURE_LIST ${CXX11_FEATURE_LIST} CACHE STRING "C++11 feature support list")
mark_as_advanced(FORCE CXX11_FEATURE_LIST)
