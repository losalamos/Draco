dnl-------------------------------------------------------------------------dnl
dnl draco/autodoc configure.in
dnl Thomas M. Evans
dnl Fri Sep  3 13:04:05 1999
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for autodoc
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([mainpage.dcc.in])
AC_CONFIG_AUX_DIR(../config)

# draco version
AC_DRACO_PREREQ
AC_DRACO_ENV
package="Draco"


dnl
dnl doxygen setup
dnl

# check for and as sign the path to doxygen
AC_PATH_PROG(DOXYGEN_PATH, doxygen, null)   

#
# Run the generic package-level autodoc configuration:
#
AC_PACKAGE_AUTODOC

#
# Do Draco-specific configuration:
#

# dotpath (for draco levelization diagram)
dotpath="${srcdir}/../doc/design"
AC_SUBST(dotpath)

# find the release number
number='null'
AC_MSG_CHECKING("draco release number")
number=`/bin/sh ../tools/draco --version`
AC_MSG_RESULT($number)
AC_SUBST(number)

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/autodoc configure.ac
dnl-------------------------------------------------------------------------dnl
