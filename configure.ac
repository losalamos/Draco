dnl-------------------------------------------------------------------------dnl
dnl draco/ configure.in
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for draco/
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([README.draco])
AC_CONFIG_AUX_DIR(config)

# draco version
AC_DRACO_PREREQ

dnl The only reason for vendor setup at this level is to produce the
dnl correct message (configure --help).   With out including this
dnl macro, none of the --with-<vendor> messages will show up.
AC_ALL_VENDORS_SETUP
AC_DRACO_ARGS

# find common tools
AC_DRACO_CHECK_TOOLS

# define fortran
AC_WITH_F90

# The order of packages in AC_CONFIG_SUBDIRS **is** important.
# The packages may depend on other packages to the left of them
# in the list.

AC_CONFIG_SUBDIRS(tools src autodoc)

# make the default head directory makefile
AC_HEAD_MAKEFILE

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/ configure.in
dnl-------------------------------------------------------------------------dnl
