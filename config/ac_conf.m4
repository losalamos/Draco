dnl-------------------------------------------------------------------------dnl
dnl ac_conf.m4
dnl service macros used in configure.in's throughout DRACO
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS
dnl
dnl add DRACO-dependent libraries necessary for a package
dnl usage: configure.in
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS, [dnl
   if test ${has_libdir:=no} != "yes" ; then
       DRACO_LIBS="${DRACO_LIBS} -L\${libdir}"
       has_libdir="yes"
   fi

   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_depends="\${libdir}/lib${lib}\${libsuffix}"
       DRACO_DEPENDS="${DRACO_DEPENDS} ${draco_depends}"
       DRACO_LIBS="${DRACO_LIBS} -l${lib}"
   done
])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_TEST
dnl
dnl add DRACO-dependent libraries necessary for a package test
dnl usage: configure.in
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS_TEST, [dnl
   DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -L\${libdir}"
   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_test_depends="\${libdir}/lib${lib}\${libsuffix}"
       DRACO_TEST_DEPENDS="${DRACO_TEST_DEPENDS} ${draco_test_depends}"
       DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -l${lib}"
   done
])

dnl-------------------------------------------------------------------------dnl
dnl AC_RUNTESTS
dnl
dnl add DRACO-package tests (default to use DejaGnu)
dnl usage: in configure.in:
dnl AC_RUNTESTS(testexec1 testexec2 ... , {nprocs1 nprocs2 ... | scalar})
dnl where serial means run as serial test only.
dnl If compiling with scalar c4 then nprocs are ignored.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_RUNTESTS, [dnl
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

AC_DEFUN(AC_TESTEXE, [dnl
   test_exe="$1"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_EXECUTABLE
dnl
dnl where executables will be installed
dnl usage: configure.in
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_EXECUTABLE, [ dnl
   install_executable="\${bindir}/\${package}"
   installdirs="${installdirs} \${bindir}"
   alltarget="${alltarget} bin/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_LIB
dnl
dnl where libraries will be installed
dnl usage: configure.in
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_LIB, [ dnl
   install_lib="\${libdir}/lib\${package}\${libsuffix}"
   installdirs="${installdirs} \${libdir}"
   alltarget="${alltarget} lib\${package}\${libsuffix}"

   # test will need to link this library
   PKG_DEPENDS='../lib${package}${libsuffix}'
   PKG_LIBS='-L.. -l${package}'
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_HEADERS
dnl
dnl where headers will be installed 
dnl usage: configure.in
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_HEADERS, [ dnl
   install_headers="\${installheaders}"
   installdirs="${installdirs} \${includedir} \${includedir}/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_BUILD_CONFIGURES
dnl
dnl build (or check if they exist) configure scripts in the current 
dnl directories subdirectories
dnl 
dnl the argument is the relative location of the draco/config
dnl directory, ie, in draco/src/configure.in the argument would be
dnl AC_BUILD_CONFIGURE(..)
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_BUILD_CONFIGURES, [ dnl

# check for configure
AC_CHECK_PROG(AUTOCONF, autoconf, autoconf, null)
if test "${AUTOCONF}" = null ; then
   AC_MSG_ERROR("Cannot find autoconf; cannot build configure scripts!")
fi

CURRENT_DIR=`pwd`

# build configure scripts in subdirs
for dir in $subdirs
   do
       if test -d ${srcdir}/${dir} ; then
	   AC_MSG_CHECKING("configure in ${srcdir}/${dir}")
	   cd $srcdir/$dir
	   if test -f configure ; then 
	       AC_MSG_RESULT("exists")
	   elif test -f configure.in ; then
	       ${AUTOCONF} --localdir=$1/../config
	       AC_MSG_RESULT("built")
	   else
	       AC_MSG_RESULT("no configure.in - skipped")
	   fi
	   cd $CURRENT_DIR
       fi
   done

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DEFINE_EIGHT_BYTE_INT_TYPE
dnl
dnl SET THE VALUE OF def_eight_byte_int_type to true, this tells
dnl dracoenv to define the EIGHT_BYTE_INT_TYPE macro to the appropriate C++
dnl data type
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DEFINE_EIGHT_BYTE_INT_TYPE, [dnl

   def_eight_byte_int_type='true'

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DEFINE_FOUR_BYTE_INT_TYPE
dnl
dnl SET THE VALUE OF def_eight_byte_int_type to true, this tells
dnl dracoenv to define the FOUR_BYTE_INT_TYPE macro to the appropriate C++
dnl data type
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DEFINE_FOUR_BYTE_INT_TYPE, [dnl

   def_four_byte_int_type='true'

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DEFINE_FOUR_BYTE_FLOAT_TYPE
dnl
dnl SET THE VALUE OF def_eight_byte_float_type to true, this tells
dnl dracoenv to define the FOUR_BYTE_FLOAT_TYPE macro to the appropriate C++
dnl data type
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DEFINE_FOUR_BYTE_FLOAT_TYPE, [dnl

   def_four_byte_float_type='true'

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DEFINE_EIGHT_BYTE_FLOAT_TYPE
dnl
dnl SET THE VALUE OF def_eight_byte_float_type to true, this tells
dnl dracoenv to define the EIGHT_BYTE_FLOAT_TYPE macro to the appropriate C++
dnl data type
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DEFINE_EIGHT_BYTE_FLOAT_TYPE, [dnl

   def_eight_byte_float_type='true'

])

dnl-------------------------------------------------------------------------dnl
dnl AC_CHECK_TOOLS
dnl
dnl Find tools used by the build system (latex, bibtex, python, etc)
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_CHECK_TOOLS, [dnl

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

])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_conf.m4
dnl-------------------------------------------------------------------------dnl

