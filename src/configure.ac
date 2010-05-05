dnl-------------------------------------------------------------------------dnl
dnl draco/src configure.in
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for draco/src/
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([configure.ac])
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

AC_CONFIG_SUBDIRS( ds++ )

##---------------------------------------------------------------------------##
## LEVEL 2
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS( xm rng c4 units traits cdi mesh_element \
		   lapack_wrap shared_lib fpe_trap ) 
# Only build plot2D if xmgrace is available
if test -n "${vendor_grace}"; then
  AC_CONFIG_SUBDIRS( plot2D )
fi


##---------------------------------------------------------------------------##
## LEVEL 3
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS( norms timestep pcgWrap viz parser \
                   cdi_gandolf cdi_analytic meshReaders \
                   special_functions )
dnl
dnl Kelly Thompson 12 January 2005
dnl Removed cdi_eospac from list of packages to build and test because
dnl EOSPAC-5.35 is not fully compatible with Lahey 6.2.  See SF
dnl tracker at https://sourceforge-web.lanl.gov/tracker/?aid=4566.
dnl
dnl cdi_eospac \
dnl

##---------------------------------------------------------------------------##
## LEVEL 4
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS( RTT_Format_Reader quadrature )

##---------------------------------------------------------------------------##
## build the makefile in src/
##---------------------------------------------------------------------------##

AC_SRC_MAKEFILE

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/configure.in
dnl-------------------------------------------------------------------------dnl
