dnl-------------------------------------------------------------------------dnl
dnl draco/src configure.in
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for draco/src/
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Makefile.in])
AC_CONFIG_AUX_DIR(../config)

# draco version
AC_DRACO_PREREQ

AC_ALL_VENDORS_SETUP
AC_DRACO_ARGS

# The order of packages in AC_CONFIG_SUBDIRS **is** important.  Each
# AC_CONFIG_SUBDIRS corresponds to a level, level 1 at the top ->
# level 6 at the bottom.  On each level the ordering is NOT
# important. See draco/doc/design/level.ps for a visualization of the
# levelization diagram

##---------------------------------------------------------------------------##
## LEVEL 1
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(stdheaders ds++)

##---------------------------------------------------------------------------##
## LEVEL 2
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(xm rng c4 units traits cdi mesh_element quadrature \
                  plot2D lapack_wrap shared_lib)

##---------------------------------------------------------------------------##
## LEVEL 3
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(timestep pcgWrap viz cdi_eospac \
                  cdi_gandolf cdi_analytic meshReaders)

##---------------------------------------------------------------------------##
## LEVEL 4
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS(RTT_Format_Reader meshReaders_Services)

##---------------------------------------------------------------------------##
## build the src/pkg configure scripts
## (this must be done after the AC_CONFIG_SUBDIRS calls)
##---------------------------------------------------------------------------##

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/configure.in
dnl-------------------------------------------------------------------------dnl
