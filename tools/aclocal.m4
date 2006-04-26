# generated automatically by aclocal 1.7.3 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

dnl-------------------------------------------------------------------------dnl
dnl ac_dracoarg.m4
dnl
dnl Declarations of Draco configure options (with some default
dnl settings). 
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:20
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_ARGS
dnl
dnl Declaration of Draco non-vendor configure options. This macro can 
dnl be called to fill out configure help screens
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_ARGS], [dnl

   dnl
   dnl Library prefix
   dnl
     
   AC_ARG_WITH(lib-prefix,
      [  --with-lib-prefix[=library prefix]
                          give prefix to libraries (default rtt_)])

   # default for lib_prefix is rtt_
   LIB_PREFIX="${with_lib_prefix:=rtt_}"
   if test "${LIB_PREFIX}" = no ; then
       LIB_PREFIX=''
   fi

   dnl
   dnl c4 toggle (scalar by default)
   dnl

   dnl define --with-c4
   AC_ARG_WITH(c4, 
      [  --with-c4[=scalar,mpi,shmem]   
		          turn on c4 (default scalar) ])

   # give with-c4 implied argument
   if test "${with_c4:=scalar}" = yes ; then
       with_c4='scalar'
   fi

   dnl
   dnl DBC toggle
   dnl

   dnl defines --with-dbc
   AC_ARG_WITH(dbc,
      [  --with-dbc[=level]      set Design-by-Contract])
	
   if test "${with_dbc}" = yes ; then
       with_dbc='7'
   elif test "${with_dbc}" = no ; then
       with_dbc='0'
   fi
	
   dnl
   dnl SHARED versus ARCHIVE libraries
   dnl

   dnl defines --enable-shared
   AC_ARG_ENABLE(shared,
      [  --enable-shared         turn on shared libraries (.a default)])

   dnl
   dnl CHOOSE A C++ COMPILER
   dnl

   dnl defines --with-cxx
   AC_ARG_WITH(cxx,
      [  --with-cxx[=gcc,icpc,sgi,kcc,compaq,guide]                                    
                          choose a c++ compiler (defaults are machine dependent)])

   dnl the default is gcc
   if test "${with_cxx}" = yes ; then
       with_cxx='gcc'
   fi

   dnl
   dnl STATIC VERSUS DYNAMIC LINKING
   dnl

   dnl defines --enable-static-ld
   AC_ARG_ENABLE(static-ld,
      [  --enable-static-ld      use (.a) libraries if possible])

   dnl
   dnl ANSI STRICT COMPLIANCE
   dnl

   dnl defines --enable-strict-ansi
   AC_ARG_ENABLE(strict-ansi,
      [  --disable-strict-ansi   turn off strict ansi compliance])

   dnl
   dnl ONE_PER INSTANTIATION FLAG
   dnl

   dnl defines --enable-one-per
   AC_ARG_ENABLE(one-per,
      [  --disable-one-per       turn off --one_per flag])

   dnl
   dnl COMPILER OPTIMZATION LEVEL
   dnl

   dnl defines --with-opt
   AC_ARG_WITH(opt,
      [  --with-opt[=0,1,2,3]    set optimization level (0 by default)])

   if test "${with_opt}" = yes ; then
       with_opt='0'
   fi

   dnl defines --enable-debug
   AC_ARG_ENABLE(debug,
      [  --enable-debug          turn on debug (-g) option])

   dnl
   dnl POSIX SOURCE
   dnl

   dnl defines --with-posix
   AC_ARG_WITH(posix,
      [  --with-posix[=num]      give posix source (system-dependent defaults)])

   dnl
   dnl ADD TO CPPFLAGS
   dnl
   
   dnl defines --with-cppflags
   AC_ARG_WITH(cppflags,
      [  --with-cppflags[=flags] add flags to \$CPPFLAGS])

   dnl
   dnl ADD TO CXXFLAGS
   dnl
   
   dnl defines --with-cxxflags
   AC_ARG_WITH(cxxflags,
      [  --with-cxxflags[=flags] add flags to \$CXXFLAGS])

   dnl
   dnl ADD TO CFLAGS
   dnl
   
   dnl defines --with-cflags
   AC_ARG_WITH(cflags,
      [  --with-cflags[=flags]   add flags to \$CFLAGS])

   dnl
   dnl ADD TO F90FLAGS
   dnl
   
   dnl defines --with-f90flags
   AC_ARG_WITH(f90flags,
      [  --with-f90flags[=flags] add flags to \$F90FLAGS])

   dnl
   dnl ADD TO ARFLAGS
   dnl
   
   dnl defines --with-arflags
   AC_ARG_WITH(arflags,
      [  --with-arflags[=flags]  add flags to \$ARFLAGS])

   dnl
   dnl ADD TO LDFLAGS
   dnl
   
   dnl defines --with-ldflags
   AC_ARG_WITH(ldflags,
      [  --with-ldflags[=flags]  add flags to \$LDFLAGS])

   dnl 
   dnl ADD TO LIBRARIES
   dnl

   dnl defines --with-libs
   AC_ARG_WITH(libs,
      [  --with-libs=[libs]      add libs to \$LIBS])

   dnl
   dnl CHOSE BIT COMPILATION ON SGI'S
   dnl

   dnl defines --enable-32-bit
   AC_ARG_ENABLE(32-bit,
      [  --enable-32-bit         do 32-bit compilation (compiler dependent)])

   dnl defines --enable-64-bit
   AC_ARG_ENABLE(64-bit,
      [  --enable-64-bit         do 64-bit compilation (compiler dependent)])

   dnl
   dnl CHOSE MIPS INSTRUCTION SET ON SGI'S
   dnl

   dnl defines --with-mips
   AC_ARG_WITH(mips,
      [  --with-mips[=1,2,3,4]   set mips, mips4 by default (SGI ONLY)])

   if test "${with_mips}" = yes ; then
       with_mips='4'
   fi

   dnl 
   dnl Arguments for options defined in ac_instrument.m4
   dnl
   
   AC_DRACO_INSTR_ARGS

   dnl
   dnl Doxygen options
   dnl

   AC_ARG_ENABLE(latex-doc,
      [  --enable-latex-doc      build latex docs with doxygen (off by default)],
      [AC_SUBST(latex_yes_no,'YES')],
      [AC_SUBST(latex_yes_no,'NO')])

   AC_ARG_WITH(doc-output,
      [  --with-doc-output=path  build documentation in path (prefix/documentation by default)],
      [AC_SUBST(doxygen_output_top,${with_doc_output})],
      [doxygen_output_top='DEFAULT'])

   dnl end of AC_DRACO_ARGS
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoarg.m4
dnl-------------------------------------------------------------------------dnl


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
   AC_ARG_WITH(stlport,
      [  --with-stlport[[=dir]]
                          replace default STL with STLPort (off by default).
                          examines value of STLPORT_BASE_DIR.
                          Only available for g++ on Linux.])

   dnl 
   dnl Coverage Analsysis
   dnl

   dnl specify type of coverage analysis.
   AC_ARG_WITH(coverage,
      [  --with-coverage[=bullseye(default)|gcov]
                          produce coverage analysis statistics (off by default).
                          examines value of COVERAGE_BASE_DIR.
                          Only available for g++ on Linux.])

   dnl 
   dnl Memory Checkers
   dnl

   dnl specify type of memory checking to be done.
   AC_ARG_WITH(memory-check,
      [  --with-memory-check[=purify(default)|insure]
                          produce binaries that are instrumented for memory 
                          checking (off by default). examines value of
                          MEMORYCHECK_BASE_DIR.
                          Only available for g++ on Linux.])

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
       if test -d ${STLPORT_BASE_DIR:=/codes/radtran/vendors/stlport/Linux}; then
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
dnl    CXX=/usr/local/bullseye/bin/g++
dnl    CC=/usr/local/bullseye/bin/gcc
dnl    LIBS="${LIBS} -L/usr/local/bullseye/lib -lcov-noThread"
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
   
   AC_MSG_CHECKING("for Bullseye installation location")
   case $host in
   *-linux-gnu)
      if ! test -d ${COVERAGE_BASE_DIR:=/usr/local/bullseye}; then
         AC_MSG_ERROR("${COVERAGE_BASE_DIR} could not be accessed.")
      fi
      ;;
   *)
      AC_MSG_ERROR("BullseyeCoverage not supported on the ${host} platform.")
      ;;
   esac
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
   case ${CXX} in
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


dnl-------------------------------------------------------------------------dnl
dnl ac_conf.m4
dnl
dnl Service macros used in configure.ac's throughout Draco.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:19
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_PREREQ
dnl
dnl Checks the configure version
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_PREREQ], [dnl

   # we need at least autoconf 2.53 to work correctly
   AC_PREREQ(2.53)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS
dnl
dnl add DRACO-dependent libraries necessary for a package
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_NEEDS_LIBS], [dnl
   if test ${has_libdir:=no} != "yes" ; then
       DRACO_LIBS="${DRACO_LIBS} -L\${libdir}"
       has_libdir="yes"
   fi

   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_depends="\${libdir}/lib\${LIB_PREFIX}${lib}\${libsuffix}"
       DRACO_DEPENDS="${DRACO_DEPENDS} ${draco_depends}"
       DRACO_LIBS="${DRACO_LIBS} -l\${LIB_PREFIX}${lib}"
   done

   # Keep a list of component dependencies free of other tags or paths.
   DEPENDENT_COMPONENTS="$1"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_TEST
dnl
dnl add DRACO-dependent libraries necessary for a package test
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_NEEDS_LIBS_TEST], [dnl
   DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -L\${libdir}"
   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_test_depends="\${libdir}/lib\${LIB_PREFIX}${lib}\${libsuffix}"
       DRACO_TEST_DEPENDS="${DRACO_TEST_DEPENDS} ${draco_test_depends}"
       DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -l\${LIB_PREFIX}${lib}"
   done
])

dnl-------------------------------------------------------------------------dnl
dnl AC_RUNTESTS
dnl
dnl add DRACO-package tests (default to use DejaGnu)
dnl usage: in configure.ac:
dnl AC_RUNTESTS(testexec1 testexec2 ... , {nprocs1 nprocs2 ... | scalar})
dnl where serial means run as serial test only.
dnl If compiling with scalar c4 then nprocs are ignored.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_RUNTESTS], [dnl
	test_alltarget="$test_alltarget $1"
        
	test_nprocs="$2"

	if test -z "${test_nprocs}" ; then
	    AC_MSG_ERROR("No procs choosen for the tests!")
        fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_TESTEXE
dnl
dnl determines what type of executable the tests are, for example, you 
dnl can set the executable to some scripting extension, like python.
dnl the default is an executable binary
dnl options are PYTHON
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_TESTEXE], [dnl
   test_exe="$1"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_EXECUTABLE
dnl
dnl where executables will be installed
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_INSTALL_EXECUTABLE], [ dnl
   install_executable="\${bindir}/\${package}"
   installdirs="${installdirs} \${bindir}"
   alltarget="${alltarget} bin/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_LIB
dnl
dnl where libraries will be installed
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_INSTALL_LIB], [ dnl
   install_lib='${libdir}/lib${LIB_PREFIX}${package}${libsuffix}'
   installdirs="${installdirs} \${libdir}"
   alltarget="${alltarget} lib\${LIB_PREFIX}\${package}\${libsuffix}"

   # test will need to link this library
   PKG_DEPENDS='../lib${LIB_PREFIX}${package}${libsuffix}'
   PKG_LIBS='-L.. -l${LIB_PREFIX}${package}'
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_HEADERS
dnl
dnl where headers will be installed 
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_INSTALL_HEADERS], [ dnl
   install_headers="\${installheaders}"
   installdirs="${installdirs} \${includedir} \${includedir}/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_CHECK_TOOLS
dnl
dnl Find tools used by the build system (latex, bibtex, python, etc)
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_CHECK_TOOLS], [dnl

   dnl
   dnl TOOL CHECKS
   dnl
   
   dnl check for and assign the path to python
   AC_PATH_PROG(PYTHON_PATH, python, null)
   if test "${PYTHON_PATH}" = null ; then
       AC_MSG_ERROR("No valid Python found!")
   fi
   
   dnl check for and assign the path to perl
   AC_PATH_PROG(PERL_PATH, perl, null)
   if test "${PERL_PATH}" = null ; then
       AC_MSG_WARN("No valid Perl found!")
   fi

   dnl check for CVS
   AC_PATH_PROG(CVS_PATH, cvs, null)
   if test "${CVS_PATH}" = null ; then
       AC_MSG_WARN("No valid CVS found!")
   fi

   dnl check for and assign the path to ghostview
   AC_CHECK_PROGS(GHOSTVIEW, ghostview gv, null)
   if test "${GHOSTVIEW}" = null ; then
       AC_MSG_WARN("No valid ghostview found!")
   fi

   dnl check for and assign the path to latex
   AC_CHECK_PROGS(LATEX, latex, null)
   if test "${LATEX}" = null ; then
       AC_MSG_WARN("No valid latex found!")
   fi
   AC_SUBST(LATEXFLAGS)

   dnl check for and assign the path to bibtex
   AC_CHECK_PROGS(BIBTEX, bibtex, null)
   if test "${BIBTEX}" = null ; then
       AC_MSG_WARN("No valid bibtex found!")
   fi
   AC_SUBST(BIBTEXFLAGS)

   dnl check for and assign the path to xdvi
   AC_CHECK_PROGS(XDVI, xdvi, null)
   if test "${XDVI}" = null ; then
       AC_MSG_WARN("No valid xdvi found!")
   fi
   AC_SUBST(XDVIFLAGS)

   dnl check for and assign the path to xdvi
   AC_CHECK_PROGS(PS2PDF, ps2pdf, null)
   if test "${PS2PDF}" = null ; then
       AC_MSG_WARN("No valid ps2pdf found!")
   fi
   dnl AC_SUBST(PS2PDFFLAGS)

   dnl check for and assign the path to xdvi
   AC_CHECK_PROGS(DOTCMD, dot, null)
   if test "${DOTCMD}" = null ; then
       AC_MSG_WARN("No valid dot found!")
   fi
   dnl AC_SUBST(DOTCMDFLAGS)

   dnl check for and assign the path to dvips
   AC_CHECK_PROGS(DVIPS, dvips, null)
   if test "${DVIPS}" = null ; then
       AC_MSG_WARN("No valid dvips found!")
   fi
   AC_SUBST(DVIPSFLAGS)

   dnl check for and assign the path for printing (lp)
   AC_CHECK_PROGS(LP, lp lpr, null)
   if test "${LP}" = null ; then
       AC_MSG_WARN("No valid lp or lpr found!")
   fi
   AC_SUBST(LPFLAGS)

   dnl check for and assign the path for doxygen
   AC_PATH_PROG(DOXYGEN_PATH, doxygen, null)
   if test "${DOXYGEN_PATH}" = null ; then
       AC_MSG_WARN("No valid Doxygen found!")
   fi
   AC_SUBST(DOXYGEN_PATH)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ASCI_WHITE_TEST_WORK_AROUND_PREPEND
dnl
dnl changes compiler from newmpxlC to newxlC so that tests can be run
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_WHITE_TEST_WORK_AROUND_PREPEND], [dnl

   # change compiler
   if test "${CXX}" = newmpxlC; then
       white_compiler='newmpxlC'
       CXX='newxlC'
       AC_MSG_WARN("Changing to ${CXX} compiler for configure tests.")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ASCI_WHITE_TEST_WORK_AROUND_APPEND
dnl
dnl changes compiler back to newmpxlC
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_WHITE_TEST_WORK_AROUND_APPEND], [dnl

   # change compiler back
   if test "${white_compiler}" = newmpxlC; then
       CXX='newmpxlC'
       AC_MSG_WARN("Changing back to ${CXX} compiler.")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_HEAD_MAKEFILE
dnl 
dnl Builds default makefile in the head directory
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HEAD_MAKEFILE], [dnl

   AC_FIND_TOP_SRC($srcdir, package_top_srcdir)
   AC_DBS_VAR_SUBSTITUTIONS
   AC_CONFIG_FILES([Makefile:config/Makefile.head.in])

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SRC_MAKEFILE
dnl 
dnl Builds default makefile in the src directory
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SRC_MAKEFILE], [dnl

   AC_FIND_TOP_SRC($srcdir, package_top_srcdir)
   AC_DBS_VAR_SUBSTITUTIONS
   AC_CONFIG_FILES([Makefile:../config/Makefile.src.in])

])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_conf.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl ac_local.m4
dnl
dnl Macros used internally within the Draco build system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_WITH_DIR
dnl
dnl Define --with-xxx[=DIR] with defaults to an environment variable.
dnl       Usage: AC_WITH_DIR(flag, CPPtoken, DefaultValue, HelpStr)
dnl                for environment variables enter \${ENVIRONVAR} for
dnl                DefaultValue
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_WITH_DIR], [dnl

 dnl
 dnl  The following M4 macros will be expanded into the body of AC_ARG_WITH
 dnl
 dnl AC_PACKAGE is the flag with all dashes turned to underscores
 dnl AC_WITH_PACKAGE will be substituted to the autoconf shell variable
 dnl    with_xxx
 dnl AC_CMDLINE is the shell command to strip double and trailing slashes
 dnl    from directory names.

 define([AC_PACKAGE], [translit($1, [-], [_])])dnl
 define([AC_WITH_PACKAGE], [with_]AC_PACKAGE)dnl
 define([AC_CMDLINE],dnl
[echo "$]AC_WITH_PACKAGE[" | sed 's%//*%/%g' | sed 's%/$%%'])dnl

 AC_ARG_WITH($1,
   [  --with-$1[=DIR]    $4 ($3 by default)],
   if test $AC_WITH_PACKAGE != "no" ; then
      if test $AC_WITH_PACKAGE = "yes" ; then
         # following eval needed to remove possible '\' from $3
         eval AC_WITH_PACKAGE=$3
      fi

      # this command removes double slashes and any trailing slash

      AC_WITH_PACKAGE=`eval AC_CMDLINE`
      if test "$AC_WITH_PACKAGE:-null}" = "null" ; then
         { echo "configure: error: --with-$1 directory is unset" 1>&2; \
           exit 1; }
      fi
      if test ! -d $AC_WITH_PACKAGE ; then
         { echo "configure: error: $AC_WITH_PACKAGE: invalid directory" 1>&2; \
           exit 1; }
      fi

      # this sets up the shell variable, with the name of the CPPtoken,
      # and that we later will do an AC_SUBST on.
      $2="${AC_WITH_PACKAGE}/"

      # this defines the CPP macro with the directory and single slash appended.
      AC_DEFINE_UNQUOTED($2, ${AC_WITH_PACKAGE}/)dnl

      # print a message to the users (that can be turned off with --silent)

      echo "$2 has been set to $$2" 1>&6

   fi)

   AC_SUBST($2)dnl

])
	
dnl-------------------------------------------------------------------------dnl
dnl AC_VENDORLIB_SETUP(1,2)
dnl
dnl set up for VENDOR_LIBS or VENDOR_TEST_LIBS
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_VENDORLIB_SETUP], [dnl

   # $1 is the vendor_<> tag (equals pkg or test)
   # $2 are the directories added 

   if test "${$1}" = pkg ; then
       VENDOR_LIBS="${VENDOR_LIBS} $2"
   elif test "${$1}" = test ; then
       VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} $2"
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_FIND_TOP_SRC(1,2)
dnl 
dnl Find the top source directory of the package by searching upward
dnl from the argument directory. The top source directory is defined
dnl as the one with a 'config' sub-directory.
dnl
dnl Note: This function will eventually quit if the searched for
dnl directory is not above the argument. It does so when $temp_dir
dnl ceases to be a valid directory, which only seems to happen after a
dnl LOT of ..'s are added to it.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_FIND_TOP_SRC], [dnl
   
   # $1 is the component's source directory
   # $2 is the variable to store the package's main source directory in.

   temp_dir=$1
   AC_MSG_CHECKING([package top source directory])
   while test -d $temp_dir -a ! -d $temp_dir/config ; do   
       temp_dir="${temp_dir}/.."
   done
   if test -d $temp_dir; then
       $2=`cd $temp_dir; pwd;`
       AC_MSG_RESULT([$$2])
   else
       AC_MSG_ERROR('Could not find package top source directory')
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl DO VARIABLE SUBSTITUTIONS ON AC_OUTPUT
dnl
dnl These are all the variable substitutions used within the draco
dnl build system
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_VAR_SUBSTITUTIONS], [dnl

   # these variables are declared "precious", meaning that they are
   # automatically substituted, put in the configure --help, and
   # cached 
   AC_ARG_VAR(CC)dnl
   AC_ARG_VAR(CFLAGS)dnl

   AC_ARG_VAR(CXX)dnl
   AC_ARG_VAR(CXXFLAGS)dnl

   AC_ARG_VAR(LD)dnl
   AC_ARG_VAR(LDFLAGS)dnl

   AC_ARG_VAR(AR)dnl
   AC_ARG_VAR(ARFLAGS)dnl

   AC_ARG_VAR(CPPFLAGS)dnl

   # dependency rules
   AC_SUBST(DEPENDENCY_RULES)

   # other compiler substitutions
   AC_SUBST(STRICTFLAG)dnl
   AC_SUBST(PARALLEL_FLAG)dnl
   AC_SUBST(RPATH)dnl
   AC_SUBST(LIB_PREFIX)dnl

   # install program
   AC_SUBST(INSTALL)dnl
   AC_SUBST(INSTALL_DATA)dnl

   # files to install
   : ${installfiles:='${install_executable} ${install_lib} ${install_headers}'}
   AC_SUBST(installfiles)dnl
   AC_SUBST(install_executable)dnl
   AC_SUBST(install_lib)dnl
   AC_SUBST(install_headers)dnl
   AC_SUBST(installdirs)dnl

   # package libraries
   AC_SUBST(alltarget)dnl
   AC_SUBST(libsuffix)dnl
   AC_SUBST(dirstoclean)dnl
   AC_SUBST(package)dnl
   AC_SUBST(DRACO_DEPENDS)dnl
   AC_SUBST(DRACO_LIBS)dnl
   AC_SUBST(VENDOR_DEPENDS)dnl
   AC_SUBST(VENDOR_INC)dnl
   AC_SUBST(VENDOR_LIBS)dnl
   AC_SUBST(ARLIBS)dnl

   # package testing libraries
   AC_SUBST(PKG_DEPENDS)dnl
   AC_SUBST(PKG_LIBS)dnl
   AC_SUBST(DRACO_TEST_DEPENDS)dnl
   AC_SUBST(DRACO_TEST_LIBS)dnl
   AC_SUBST(VENDOR_TEST_DEPENDS)dnl
   AC_SUBST(VENDOR_TEST_LIBS)dnl
   AC_SUBST(ARTESTLIBS)dnl
   AC_SUBST(test_alltarget)dnl
   AC_SUBST(test_flags)dnl
   AC_SUBST(test_scalar)dnl
   AC_SUBST(test_nprocs)dnl
   AC_SUBST(test_output_files)dnl

   # libraries
   AC_ARG_VAR(LIBS)dnl

   # configure options
   AC_SUBST(configure_command)dnl

   # directories in source tree
   AC_SUBST(package_top_srcdir)
   
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_local.m4
dnl-------------------------------------------------------------------------dnl

