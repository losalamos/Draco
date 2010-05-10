dnl ----------------------------------------------------------------------- dnl
dnl File  : draco/config/ac_instrument.m4
dnl Author: Kelly Thompson
dnl Date  : 2006 MAR 20
dnl
dnl Defines the Draco build system environment needed for
dnl instrumentation.  Provide support for STLport and BullseyeCoverage
dnl on Linx.
dnl
dnl ----------------------------------------------------------------------- dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

dnl ----------------------------------------------------------------------- dnl
dnl AC_DRACO_INSTR_ARGS
dnl
dnl Called by : ac_dracoarg.m4
dnl Purpose   : Provide help/usage messages for the features in this file.
dnl ----------------------------------------------------------------------- dnl

AC_DEFUN([AC_DRACO_INSTR_ARGS], [dnl

   dnl 
   dnl STLport 
   dnl

   dnl Request a build that uses STLPort (specify location of STLPort).
   AC_ARG_WITH([stlport],
     [AS_HELP_STRING([--with-stlport=DIR],
       [replace default STL with STLPort (off by default). Examines
        value of @S|@STLPORT_BASE_DIR. Only available for g++ on Linux.])])

   dnl 
   dnl Coverage Analsysis
   dnl

   dnl specify type of coverage analysis.
   AC_ARG_WITH([coverage],
     [AS_HELP_STRING([--with-coverage@<:@=bullseye(default)|gcov@:>@],
       [produce coverage analysis statistics (off by
        default). Examines value of @S|@COVERAGE_BASE_DIR. Only
        available for g++ on Linux.])]) 

   dnl 
   dnl Memory Checkers
   dnl

   dnl specify type of memory checking to be done.
   AC_ARG_WITH(memory-check,
     [AS_HELP_STRING([--with-memory-check@<:@=purify(default)|insure@:>@],
       [produce binaries that are instrumented for memory checking
        (off by default). examines value of
        @S|@MEMORYCHECK_BASE_DIR. Only available for g++ on Linux.])]) 

])

dnl ----------------------------------------------------------------------- dnl
dnl AC_DRACO_INSTR_ENV
dnl
dnl Called by : ac_dracoenv.m4
dnl Purpose   : Provide a single function that can be called from 
dnl             ac_dracoarg.m4 (AC_DRACO_ENV) to modify the build 
dnl             environment if the user requests any of the instrument
dnl             options.
dnl ----------------------------------------------------------------------- dnl

AC_DEFUN([AC_DRACO_INSTR_ENV], [dnl

   # we must know the host
   AC_REQUIRE([AC_CANONICAL_HOST])

   AC_DBS_STLPORT_ENV
   AC_DBS_COVERAGE_ENV
   AC_DBS_MEMORY_CHECK_ENV

])
 
dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_STLPORT_ENV
dnl
dnl Used by AC_DRACO_ENV, this macro checks the configure line for the
dnl presence of "--with-stlport".  If this option is found, the build
dnl system's environment is modified so that all the all C++ compiles
dnl use the STL libraries included with STLPort instead of the
dnl compiler's native STL defintions.
dnl If --with-stlport is on the configure line, we must prepend
dnl CXXFLAGS and CPPFLAGS with -I<path_to_stlport>.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_STLPORT_ENV], [dnl

   AC_MSG_CHECKING("option: STLPort?")
   AC_MSG_RESULT("${with_stlport:=no}")

   # Provide an error if this is not Linux
   if test ${with_stlport} != no; then
     case ${host} in
     *-linux-gnu)
       ;;
     *)
       AC_MSG_ERROR("STLPort not supported on the ${host} platform.")
       ;;
     esac
   fi

   if test ${with_stlport} != no; then

     # Find STLPort's location
     AC_MSG_CHECKING("for STLPort installation location")

     # if --with-stlport is requested with no dir specified, then check
     # the value of STLPORT_BASE_DIR.
     if test ${with_stlport} = yes; then
       if test -d ${STLPORT_BASE_DIR:=/ccs/codes/radtran/vendors/stlport/Linux}; then
         with_stlport=${STLPORT_BASE_DIR}
       else
         AC_MSG_ERROR("${STLPORT_BASE_DIR} could not be accessed.")
       fi
     fi
     AC_MSG_RESULT("${with_stlport}")
  
     # Double check accessibility.
  
     if ! test -d "${with_stlport}/include"; then
        AC_MSG_ERROR("Invalid directory $with_stlport}/include")
     fi
     if ! test -r "${with_stlport}/lib/libstlportstlg.so"; then
        AC_MSG_ERROR("Invalid library ${with_stlport}/lib/libstlportstlg.so")
     fi
  
     # Modify environment
  
     AC_MSG_CHECKING("STLPort modification for CPPFLAGS")
     cppflag_mods="-I${with_stlport}/include -D_STLP_DEBUG"
     dnl Consider adding -D_STLP_DEBUG_UNINITIALIZED
     CPPFLAGS="${cppflag_mods} ${CPPFLAGS}"
     AC_MSG_RESULT([${cppflag_mods}])
  
  dnl Problems with STLport-5.0.X prevent us from using the optimized specializations.
  
     AC_MSG_CHECKING("STLPort modification for LIBS")
     libs_mods="-L${with_stlport}/lib -lstlportstlg"
     LIBS="${libs_mods} ${LIBS}"
     AC_MSG_RESULT([${libs_mods}])
  
     AC_MSG_CHECKING("STLPort modifications for RPATH")
     rpath_mods="-Xlinker -rpath ${with_stlport}/lib"
     RPATH="${rpath_mods} ${RPATH}"
     AC_MSG_RESULT("$rpath_mods}")

   fi dnl  if test ${with_stlport} != no

   dnl end of AC_DBS_STLPORT_ENV
])


dnl ------------------------------------------------------------------------dnl
dnl AC_DBS_COVERAGE_ENV
dnl
dnl Used by AC_DRACO_ENV, this macro checks the configure line for the
dnl presence of "--with-coverage[=<bullseye|gcov>]".  If this option
dnl is found, the build system's environment is modified so that all
dnl the all C++ compiles use the compilers provided by the coverage
dnl tool and the coverage tool's libraries must be added to the list
dnl of LIBS.
dnl
dnl If support for another coverage tool is added here, then the main
dnl body of code needs to be replaced with a case statement for each
dnl tool.  The environment modification for each tool should be in its
dnl own function.
dnl
dnl Defines:
dnl    with_coverage
dnl    COVERAGE_BASE_DIR
dnl
dnl Modifies:
dnl    CXX, CC, LIBS
dnl
dnl Bullseye specifics:
dnl
dnl If --with-coverage[=bulleye] is on the configure line, we must set:
dnl    CXX=/ccs/opt/x86/bullseye/bin/g++
dnl    CC=/ccs/opt/x86/bullseye/bin/gcc
dnl    LIBS="${LIBS} -L/ccs/opt/x86/bullseye/lib -lcov-noThread"
dnl ------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_COVERAGE_ENV], [dnl

   AC_MSG_CHECKING("Option: coverage analysis?")
   if test "${with_coverage:=no}" = no; then
      AC_MSG_RESULT("${with_coverage}")      
   else
      case ${with_coverage} in
      [bB]ullseye | BULLSEYE | yes )
         with_coverage=bullseye
         AC_MSG_RESULT("${with_coverage}")
         AC_DBS_BULLSEYE_ENV
      ;;
      gcov)
         AC_MSG_ERROR("Support for gcov has not been implemented.")
      ;;
      *)
         AC_MSG_ERROR("Unknown coverage tool ${with_coverage}.")
      ;;
      esac
   fi

   dnl end of AC_DBS_COVERAGE_ENV
])

dnl ------------------------------------------------------------------------dnl
dnl AC_DBS_BULLSEYE_ENV
dnl
dnl Modify build environment to support BullseyeCoverage analsysis (Linux
dnl only). 
dnl ------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_BULLSEYE_ENV], [dnl

   # Check availability
   
   case $host in
   *-linux-gnu)
      AC_PATH_PROG( BULLSEYECOVERAGE, cov01, null )
      COVERAGE_BASE_DIR=`echo ${BULLSEYECOVERAGE} | sed -e 's/\/bin\/.*//'`
      if ! test -d ${COVERAGE_BASE_DIR}; then
         AC_MSG_ERROR("${COVERAGE_BASE_DIR} could not be accessed.")
      fi
      ;;
   *)
      AC_MSG_ERROR("BullseyeCoverage not supported on the ${host} platform.")
      ;;
   esac
   AC_MSG_CHECKING("for Bullseye installation location")
   AC_MSG_RESULT("${COVERAGE_BASE_DIR}")

   # Double check accessibility and other requirements

   if ! test -d "${COVERAGE_BASE_DIR}/include"; then
      AC_MSG_ERROR("Invalid directory ${COVERAGE_BASE_DIR}/include")
   fi
   if ! test -r "${COVERAGE_BASE_DIR}/lib/libcov-noThread.a"; then
      AC_MSG_ERROR("Invalid library ${COVERAGE_BASE_DIR}/lib/libcov-noThread.a")
   fi
   if ! test -x "${COVERAGE_BASE_DIR}/bin/cov01"; then
      AC_MSG_ERROR("Couldn't execute ${COVERAGE_BASE_DIR}/bin/cov01")
   fi

   # BullseyeCoverage only works with g++, gcc, icc, and icpc

   AC_MSG_CHECKING("Bullseye equivalent compiler")
   short_cxx=`echo ${CXX} | sed -e 's/.*\///g'`
   case ${short_cxx} in
   g++ | gcc)
      CXX=${COVERAGE_BASE_DIR}/bin/g++
      CC=${COVERAGE_BASE_DIR}/bin/gcc
      AC_MSG_RESULT("${CXX}")
      ;;
   icc | icpc)
      CXX=${COVERAGE_BASE_DIR}/bin/icpc
      CC=${COVERAGE_BASE_DIR}/bin/icc
      AC_MSG_RESULT("${CXX}")
      ;;
   *)
      AC_MSG_ERROR("CXX must be one of g++, gcc, icc or icpc")
      ;;
   esac

   # Modify environment

   AC_MSG_CHECKING("Bullseye modification for LIBS")
   libs_mods="-L${COVERAGE_BASE_DIR}/lib -lcov-noThread"
   LIBS="${libs_mods} ${LIBS}"
   AC_MSG_RESULT([${libs_mods}])

# Turn of DBC checks at these screw up coverage numbers.
   if test "${with_dbc:-yes}" != 0; then
      with_dbc=0
      AC_MSG_WARN("Design-by-Contract assertions have been disabled for due to activation of code coverage mode.")
   fi

   dnl end of AC_DBS_COVERAGE_ENV
])

dnl ------------------------------------------------------------------------dnl
dnl AC_DBC_MEMORY_CHECK_ENV
dnl
dnl Modify environemnt to support memory profiling via Purify, Insure++, etc.
dnl ------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_MEMORY_CHECK_ENV], [dnl

  AC_MSG_CHECKING("Option: memory checking?")
  if test "${with_memory_check:=no}" = no; then
     AC_MSG_RESULT("none")
  else
     AC_MSG_ERROR("This feature is not enabled at this time.")
  fi

   dnl end of AC_DBS_MEMORY_CHECK_ENV
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_instrument.m4
dnl-------------------------------------------------------------------------dnl

