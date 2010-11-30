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

dnl VENDORS:
dnl - This call is needed here so that the local configure script
dnl   knows what vendors are available for the current built.
dnl - The design works closely with the 'modules' feature provided by
dnl   draco/environment/Modules.  For example, loading the openmpi
dnl   module will set the environment variables MPI_INC_DIR and
dnl   MPI_LIB_DIR.   The configure script looks for these environment
dnl   variables to automatically configure the MPI vendor and package
dnl   c4. 
dnl - Alternatively, you can specify vendor configurations on the
dnl   command line using --with-<vendor>=library|keyword,
dnl   --with-<vendor>-inc=<DIR>, and --with-<vendor>-lib=<DIR>.  For
dnl   example, --with-mpi=openmpi
dnl   --with-mpi-inc=/ccs/codes/mpi/openmpi/Linux64/1.3.3/include
dnl   --with-mpi-lib=/ccs/codes/mpi/openmpi/Linux64/1.3.3/lib
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

dnl Warning: if we do not find an MPI vendor, c4 will default to
dnl scalar-mode.
AC_MPI_SETUP(pkg)
AC_CONFIG_SUBDIRS( xm rng c4 units traits cdi mesh_element \
 		   shared_lib fpe_trap ode linear ) 

dnl Determine if these vendors are available.
AC_LAPACK_SETUP(pkg)
AC_GRACE_SETUP(pkg)
AC_GSL_SETUP(pkg)

# Only build lapack_wrap if lapack is available
if test -d ${srcdir}/lapack_wrap; then
  if test -n "${vendor_lapack}"; then
    AC_CONFIG_SUBDIRS( lapack_wrap )
  else
    AC_NO_VENDOR_WARN([lapack_wrap],[lapack])
  fi
fi

# Only build plot2D if xmgrace is available
if test -d ${srcdir}/plot2D; then
  if test -n "${vendor_grace}"; then
    AC_CONFIG_SUBDIRS( plot2D )
  else
    AC_NO_VENDOR_WARN([plot2D],[grace])
  fi
fi

##---------------------------------------------------------------------------##
## LEVEL 3
##---------------------------------------------------------------------------##

AC_CONFIG_SUBDIRS( norms timestep diagnostics viz parser \
   cdi_analytic meshReaders roots fit )

AC_GANDOLF_SETUP(pkg)

if test -d ${srcdir}/cdi_gandolf; then
  if test -n "${vendor_gandolf}"; then
    AC_CONFIG_SUBDIRS( cdi_gandolf )
  else
    AC_NO_VENDOR_WARN([cdi_gandolf],[gandolf])
  fi
fi

# Only build pcgWrap if pcg and lapack are available
if test -d ${srcdir}/pcgWrap && ! test "${with_pcg:-no}" = "no"; then
  if test -n "${vendor_lapack}" && test -n "${vendor_pcg}"; then
    AC_CONFIG_SUBDIRS( pcgWrap )
  else
    AC_NO_VENDOR_WARN([pcgWrap],[pcg])
  fi
fi

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

AC_CONFIG_SUBDIRS( RTT_Format_Reader min )

# Only build special_functions if GSL is available.
if test -d ${srcdir}/special_functions; then
  if test -n "${vendor_gsl}"; then
    AC_CONFIG_SUBDIRS( special_functions )
  else
    AC_NO_VENDOR_WARN([special_functions],[gsl])
  fi
fi

##---------------------------------------------------------------------------##
## LEVEL 5
##---------------------------------------------------------------------------##

# Only build quadrature if GSL is available.
if test -d ${srcdir}/quadrature; then
  if test -n "${vendor_gsl}"; then
    AC_CONFIG_SUBDIRS( quadrature )
  else
    AC_NO_VENDOR_WARN([quadrature],[gsl])
  fi
fi

##---------------------------------------------------------------------------##
## build the makefile in src/
##---------------------------------------------------------------------------##

AC_SRC_MAKEFILE

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/configure.in
dnl-------------------------------------------------------------------------dnl
