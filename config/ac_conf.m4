dnl-------------------------------------------------------------------------dnl
dnl ac_conf.m4
dnl
dnl Service macros used in configure.ac's throughout Draco.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:19
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

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

   newtests="$1"
   test_nprocs="$2"

   if test -z "${test_nprocs}" ; then
     AC_MSG_ERROR("No procs choosen for the tests!")
   fi

dnl If we ran AC_RUNTESTS with "serial" then mark it so here.
     if test "$test_nprocs" = "serial" || test "$test_nprocs" = "scalar" ; then
       scalar_tests="$scalar_tests $newtests"   
     else
       parallel_tests="$parallel_tests $newtests"
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
dnl AC_TEST_APPLICATION
dnl
dnl Register a binary to be tested during "make test".  This function
dnl takes two arguments. The first argument is a space delimited list
dnl of binaries to be tested.  These second argument is the number of
dnl processors to use in these tests.  The 
dnl
dnl All of tests are run as parallel tests.
dnl
dnl Example use in configure.ac:
dnl
dnl AC_TEST_APPLICATION(tstSerranoExe \
dnl                     tstJibberish \
dnl                     , 1 2 6)
dnl
dnl Background:
dnl

dnl The new m4 macro AC_TEST_APPLICATION has been implemented to
dnl support a specific type of unit test on BPROC systems.  The
dnl Capsaicin project has unit tests written in C++ that create new
dnl processes that run under MPI.  For example, the unit test
dnl tstSerranoExe will setup an input deck and execute 
dnl "mpirun -np 4 ../bin/serrano test1.inp". Output is captured and
dnl evaluated by tstSerranoExe. On most systems, this type of unit
dnl test could be executed as a normal unit test.  However, BProc
dnl systems do not allow mpirun to be executed from the backend. Thus,
dnl AC_TEST_APPLICATION tests will execute the unit test as a scalar
dnl test with the expectation that an mpi process will be started by
dnl the unit test.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_TEST_APPLICATION], [dnl

   if test -z "$2" ; then
     AC_MSG_ERROR("ac_test_application requires 2 arguments.")
   fi

   if test ! "${app_test_nprocs:-none}" = "none" ; then
     AC_MSG_WARN("More than one call to ac_test_application, using nproc info from last call!")
   fi   

   app_test_nprocs="$2"
   app_tests="$app_tests $1"

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
dnl AC_ASCI_PURPLE_TEST_WORK_AROUND_PREPEND
dnl
dnl changes compiler from newmpxlC to newxlC so that tests can be run
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_PURPLE_TEST_WORK_AROUND_PREPEND], [dnl

   # change compiler
   case ${CXX} in
   */newmpxlC)
       purple_compiler='newmpxlC'
       CXX='newxlC'
       AC_MSG_WARN("Changing to ${CXX} compiler for configure tests.")
       ;;
   esac

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ASCI_PURPLE_TEST_WORK_AROUND_APPEND
dnl
dnl changes compiler back to newmpxlC
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_PURPLE_TEST_WORK_AROUND_APPEND], [dnl

   # change compiler back
   if test "${purple_compiler}" = newmpxlC; then
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

