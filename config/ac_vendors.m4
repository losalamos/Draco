dnl-------------------------------------------------------------------------dnl
dnl ac_vendors.m4
dnl
dnl Macros for each vendor that is used supported by the Draco build
dnl system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## All vendor macros should take the following arguments:
##     pkg      - this vendor is used in the package (default)
##     test     - this vendor is used only in the package's test
##
## Each vendor requires an AC_<VENDOR>_SETUP and AC_<VENDOR>_FINALIZE
## function.
##---------------------------------------------------------------------------##

dnl-------------------------------------------------------------------------dnl
dnl AC_MPI_SETUP
dnl
dnl MPI implementation (off by default)
dnl MPI is an optional vendor
dnl
dnl we wait to set the basic MPI libraries (if it is on) until
dnl after checking the C4 status; these functions are performed
dnl in ac_dracoenv.m4, section SYSTEM-SPECIFIC SETUP; we do this
dnl here because each platform has different mpi options for
dnl vendors and mpich
dnl
dnl note that we used to do this in a function called AC_COMM_SET;
dnl however, there are too many platform-dependent variables 
dnl to continue doing this; so we do all these operations in the
dnl platform specific section of ac_dracoenv.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_MPI_SETUP], [dnl

   dnl define --with-mpi
   AC_ARG_WITH(mpi,
      [  --with-mpi=[vendor,mpich] 
	                  determine MPI implementation (vendor on SGI,SUN; mpich on LINUX)])

   dnl define --with-mpi-inc and --with-mpi-lib
   AC_WITH_DIR(mpi-inc, MPI_INC, \${MPI_INC_DIR},
	       [tell where MPI includes are])
   AC_WITH_DIR(mpi-lib, MPI_LIB, \${MPI_LIB_DIR},
	       [tell where MPI libs are])

   # determine if this package is needed for testing or for the
   # package
   vendor_mpi=$1 

   # set default value for with_mpi which is no
   if test "${with_mpi:=no}" = yes ; then 
       with_mpi='vendor'
   fi

   # if the user sets MPI_INC and MPI_LIB directories then turn on  
   # with_mpi and set it to vendor if with_mpi=no to begin with
   if test "${with_mpi}" = no ; then
       if test -n "${MPI_INC}" ; then
	   with_mpi='vendor'
       elif test -n "${MPI_LIB}" ; then
	   with_mpi='vendor'
       fi
   fi

]) 

##---------------------------------------------------------------------------##

AC_DEFUN([AC_MPI_FINALIZE], [dnl

   # only add stuff if mpi is not no and the vendor is defined
   if test "${with_mpi}" != no && test -n "${vendor_mpi}"; then

       # include path
       if test -n "${MPI_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${MPI_INC}"
       fi
   
       # libraries
       if test -n "${MPI_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} ${mpi_libs})
       elif test -z "${MPI_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_mpi, ${mpi_libs})
       fi

       # add MPI directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${MPI_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${MPI_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPRNG_SETUP
dnl
dnl SPRNG LIBRARY SETUP (on by default -lfg)
dnl SPRNG is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SPRNG_SETUP], [dnl

   dnl define --with-sprng
   AC_ARG_WITH(sprng,
      [  --with-sprng[=lib]      determine the rng lib (lfg is default)])
	
   dnl define --with-sprng-inc and --with-sprng-lib
   AC_WITH_DIR(sprng-inc, SPRNG_INC, \${SPRNG_INC_DIR},
	       [tell where SPRNG includes are])
   AC_WITH_DIR(sprng-lib, SPRNG_LIB, \${SPRNG_LIB_DIR},
	       [tell where SPRNG libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_sprng=$1

   # choices are with_sprng = lfg, lcg, yes, or no

   # default (sprng is set to lfg by default)
   if test "${with_sprng:=lfg}" = yes ; then
       with_sprng='lfg'
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_SPRNG_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_sprng}"; then

       # include path
       if test -n "${SPRNG_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${SPRNG_INC}"
       fi
   
       # libraries
       if test -n "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -L${SPRNG_LIB} -l${with_sprng})
       elif test -z "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -l${with_sprng})
       fi

       # add sprng directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SPRNG_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SPRNG_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_AZTEC_SETUP
dnl
dnl AZTEC SETUP (on by default)
dnl AZTEC is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_AZTEC_SETUP], [dnl

   dnl define --with-aztec
   AC_ARG_WITH(aztec,
      [  --with-aztec=[lib]      determine the aztec lib (aztec is the default])
 
   dnl define --with-aztec-inc
   AC_WITH_DIR(aztec-inc, AZTEC_INC, \${AZTEC_INC_DIR},
	       [tell where AZTEC includes are])

   dnl define --with-aztec-lib
   AC_WITH_DIR(aztec-lib, AZTEC_LIB, \${AZTEC_LIB_DIR},
	       [tell where AZTEC libraries are])

   # set default value of aztec includes and libs
   if test "${with_aztec:=aztec}" = yes ; then
       with_aztec='aztec'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_aztec=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_AZTEC_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_aztec}" ; then

       # include path
       if test -n "${AZTEC_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${AZTEC_INC}"
       fi

       # library path
       if test -n "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -L${AZTEC_LIB} -l${with_aztec})
       elif test -z "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -l${with_aztec})
       fi

       # add AZTEC directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${AZTEC_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${AZTEC_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GSL_SETUP
dnl
dnl GSL SETUP (on by default)
dnl GSL is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_GSL_SETUP, [dnl

   dnl define --with-gsl
   AC_ARG_WITH(gsl,
      [  --with-gsl=[lib]      determine the gsl lib (gsl is the default])
 
   dnl define --with-gsl-inc
   AC_WITH_DIR(gsl-inc, GSL_INC, \${GSL_INC_DIR},
	       [tell where GSL includes are])

   dnl define --with-gsl-lib
   AC_WITH_DIR(gsl-lib, GSL_LIB, \${GSL_LIB_DIR},
	       [tell where GSL libraries are])

   # set default value of gsl includes and libs
   if test "${with_gsl:=gsl}" = yes ; then
       with_gsl='gsl'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_gsl=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_GSL_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_gsl}"; then

       # include path
       if test -n "${GSL_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GSL_INC}"
       fi

       # library path
       if test -n "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, -L${GSL_LIB} -l${with_gsl})
       elif test -z "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, -l${with_gsl})
       fi

       # add GSL directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GSL_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GSL_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GSLCBLAS_SETUP
dnl
dnl GSLCBLAS SETUP (on by default)
dnl GSLCBLAS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GSLCBLAS_SETUP], [dnl

   dnl define --with-gslcblas
   AC_ARG_WITH(gslcblas,
      [  --with-gslcblas=[lib]      determine the gslcblas lib (gslcblas is the default])
 
   dnl define --with-gslcblas-inc
   AC_WITH_DIR(gslcblas-inc, GSLCBLAS_INC, \${GSLCBLAS_INC_DIR},
	       [tell where GSLCBLAS includes are])

   dnl define --with-gslcblas-lib
   AC_WITH_DIR(gslcblas-lib, GSLCBLAS_LIB, \${GSLCBLAS_LIB_DIR},
	       [tell where GSLCBLAS libraries are])

   # set default value of gslcblas includes and libs
   if test "${with_gslcblas:=gslcblas}" = yes ; then
       with_gslcblas='gslcblas'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_gslcblas=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_GSLCBLAS_FINALIZE], [dnl

   # set up the libraries and include path
   if test "${vendor_gslcblas}"; then

       # include path
       if test -n "${GSLCBLAS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GSLCBLAS_INC}"
       fi

       # library path
       if test -n "${GSLCBLAS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gslcblas, -L${GSLCBLAS_LIB} -l${with_gslcblas})
       elif test -z "${GSLCBLAS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gslcblas, -l${with_gslcblas})
       fi

       # add GSLCBLAS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GSLCBLAS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GSLCBLAS_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_TRILINOS_SETUP
dnl
dnl TRILINOS SETUP (on by default)
dnl TRILINOS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_TRILINOS_SETUP], [dnl

   dnl define --with-trilinos
   AC_ARG_WITH(trilinos,
      [  --with-trilinos=[lib]    determine the trilinos implementation (aztecoo is default])
 
   dnl define --with-trilinos-inc
   AC_WITH_DIR(trilinos-inc, TRILINOS_INC, \${TRILINOS_INC_DIR},
	       [tell where TRILINOS includes are])

   dnl define --with-trilinos-lib
   AC_WITH_DIR(trilinos-lib, TRILINOS_LIB, \${TRILINOS_LIB_DIR},
	       [tell where TRILINOS libraries are])

   # set default value of trilinos includes and libs
   if test "${with_trilinos:=aztecoo}" = yes ; then
       with_trilinos='aztecoo'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_trilinos=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_TRILINOS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_trilinos}" ; then

       # include path
       if test -n "${TRILINOS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${TRILINOS_INC}"
       fi

       # library path
       if test -n "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -L${TRILINOS_LIB} -l${with_trilinos} -lepetra -ltriutils -ly12m)
       elif test -z "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -l${with_trilinos} -lepetra -ltriutils -ly12m)
       fi

       # add TRILINOS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${TRILINOS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${TRILINOS_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_PCG_SETUP
dnl
dnl PCG LIBRARY SETUP (on by default)
dnl PCG is a required vendor
dnl
dnl note that we add some system-specific libraries for this
dnl vendor in AC_DRACO_ENV; also, the user must set up LAPACK for
dnl this to work
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_PCG_SETUP], [dnl

   dnl define --with-pcg
   AC_ARG_WITH(pcg,        
      [  --with-pcg[=lib]        determine the pcg lib name (pcg is default)])

   dnl define --with-pcg-lib
   AC_WITH_DIR(pcg-lib, PCG_LIB, \${PCG_LIB_DIR},
	       [tell where PCG libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_pcg=$1

   # pcg is set to libpcg by default
   if test "${with_pcg:=pcg}" = yes ; then
       with_pcg='pcg'
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_PCG_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_pcg}"; then

       # library path
       if test -z "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -l${with_pcg})
       elif test -n "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -L${PCG_LIB} -l${with_pcg})
       fi

       # add PCG directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${PCG_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GANDOLF_SETUP
dnl
dnl GANDOLF LIBRARY SETUP (on by default)
dnl GANDOLF is a required vendor
dnl
dnl SGI needs "-lfortran" on the load line when including libgandolf.a.
dnl This library is added to ${LIBS} in AC_DRACO_ENV.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GANDOLF_SETUP], [dnl

   dnl define --with-gandolf
   AC_ARG_WITH(gandolf,        
      [  --with-gandolf[=lib]    determine the gandolf lib name (gandolf is default)])

   dnl define --with-gandolf-lib
   AC_WITH_DIR(gandolf-lib, GANDOLF_LIB, \${GANDOLF_LIB_DIR},
	       [tell where GANDOLF libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_gandolf=$1

   # gandolf is set to libgandolf by default
   if test "${with_gandolf:=gandolf}" = yes ; then
       with_gandolf='gandolf'
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_GANDOLF_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_gandolf}"; then

       # set up library paths
       if test -z "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -l${with_gandolf})
       elif test -n "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -L${GANDOLF_LIB} -l${with_gandolf})
       fi

       # add GANDOLF directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GANDOLF_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_EOSPAC5_SETUP
dnl
dnl EOSPAC5 LIBRARY SETUP (on by default)
dnl EOSPAC5 is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_EOSPAC5_SETUP], [dnl

   dnl define --with-eospac
   AC_ARG_WITH(eospac,        
      [  --with-eospac[=lib]     determine the eospac lib name (eospac is default)])

   dnl define --with-eospac-lib
   AC_WITH_DIR(eospac-lib, EOSPAC5_LIB, \${EOSPAC5_LIB_DIR},
	       [tell where EOSPAC5 libraries are])

   # determine if this package is needed for testing or for the 
   # package (valid values are pkg or test)
   vendor_eospac=$1

   # eospac is set to libeospac by default
   if test "${with_eospac:=eospac}" = yes ; then
       with_eospac='eospac'
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_EOSPAC5_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_eospac}"; then

       # set up library paths
       if test -z "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -l${with_eospac})
       elif test -n "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -L${EOSPAC5_LIB} -l${with_eospac})
       fi

       # add EOSPAC5 directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${EOSPAC5_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_LAPACK_SETUP
dnl
dnl LAPACK SETUP (on by default)
dnl LAPACK is a required vendor
dnl
dnl NOTE: this also sets up the BLAS
dnl
dnl note that we add system specific libraries to this list in
dnl ac_platforms.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_LAPACK_SETUP], [dnl

   dnl define --with-lapack
   AC_ARG_WITH(lapack,
      [  --with-lapack=[vendor,atlas]
                          determine LAPACK implementation (vendor default)])

   dnl define --with-lapack-lib
   AC_WITH_DIR(lapack-lib, LAPACK_LIB, \${LAPACK_LIB_DIR}, 
	       [tell where LAPACK libs are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_lapack=$1

   # lapack is set to vendor by default
   if test "${with_lapack:=vendor}" = yes ; then
       with_lapack='vendor'
   fi

   # define the atlas libraries (these are system independent)
   if test "${with_lapack}" = atlas; then
       lapack_libs='-llapack -lf77blas -lcblas -latlas'
   fi
])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_LAPACK_FINALIZE], [dnl

   # set up lapack libraries
   if test -n "${vendor_lapack}"; then

       # set libraries
       if test -z "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, ${lapack_libs})
       elif test -n "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} ${lapack_libs})
       fi

       # add LAPACK directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${LAPACK_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GRACE_SETUP
dnl
dnl GRACE SETUP (on by default)
dnl GRACE is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GRACE_SETUP], [dnl

   dnl define --with-grace
   AC_ARG_WITH(grace,
      [  --with-grace=[lib]      determine the grace lib (grace_np is the default])
 
   dnl define --with-grace-inc
   AC_WITH_DIR(grace-inc, GRACE_INC, \${GRACE_INC_DIR},
	       [tell where GRACE includes are])

   dnl define --with-grace-lib
   AC_WITH_DIR(grace-lib, GRACE_LIB, \${GRACE_LIB_DIR},
	       [tell where GRACE libraries are])

   # set default value of grace includes and libs
   if test "${with_grace:=grace_np}" = yes ; then
       with_grace='grace_np'
   fi

   # define GRACE header file
   GRACE_H="<${with_grace}.h>"
   AC_DEFINE_UNQUOTED(GRACE_H, ${GRACE_H})dnl

   # determine if this package is needed for testing or for the 
   # package
   vendor_grace=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_GRACE_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_grace}" ; then

       # include path
       if test -n "${GRACE_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GRACE_INC}"
       fi

       # library path
       if test -n "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -L${GRACE_LIB} -l${with_grace})
       elif test -z "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -l${with_grace})
       fi

       # add GRACE directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GRACE_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GRACE_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_VENDOR_FINALIZE
dnl
dnl Run at the end of the environment setup to add defines required by
dnl the vendors.  We do this to allow platform specific mods to the 
dnl vendor defines BEFORE they are added to CCPFLAGS, etc. 
dnl
dnl This macro needs to be updated when new vendors are added.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_VENDOR_FINALIZE], [dnl

   # call finalize functions for each vendor, the order is important
   # each vendor setup is appended to the previous; thus, the calling
   # level goes from high to low
   AC_TRILINOS_FINALIZE

   AC_AZTEC_FINALIZE
   AC_PCG_FINALIZE

   AC_LAPACK_FINALIZE
   AC_EOSPAC5_FINALIZE
   AC_GANDOLF_FINALIZE
   AC_SPRNG_FINALIZE
   AC_GRACE_FINALIZE

   AC_GSL_FINALIZE
   AC_GSLCBLAS_FINALIZE

   AC_MPI_FINALIZE

   # print out vendor include paths
   AC_MSG_CHECKING("vendor include paths")
   if test -n "${VENDOR_INC_DIRS}"; then
       AC_MSG_RESULT("${VENDOR_INC_DIRS}")
   else
       AC_MSG_RESULT("no vendor include dirs defined")
   fi

   # print out vendor lib paths
   AC_MSG_CHECKING("vendor lib paths")
   if test -n "${VENDOR_LIB_DIRS}"; then
       AC_MSG_RESULT("${VENDOR_LIB_DIRS}")
   else
       AC_MSG_RESULT("no vendor lib dirs defined")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ALL_VENDORS_SETUP
dnl
dnl DRACO INCLUSIVE VENDOR MACRO
dnl-------------------------------------------------------------------------dnl
dnl allows one to include all vendor macros by calling this macro.
dnl designed for draco/configure.in and draco/src/configure.in

AC_DEFUN(AC_ALL_VENDORS_SETUP, [dnl

   dnl include all macros for easy use in top-level configure.in's
   AC_MPI_SETUP(pkg)
   AC_SPRNG_SETUP(pkg)
   AC_PCG_SETUP(pkg)
   AC_AZTEC_SETUP(pkg)
   AC_GSL_SETUP(pkg)
   AC_GSLCBLAS_SETUP(pkg)
   AC_TRILINOS_SETUP(pkg)
   AC_LAPACK_SETUP(pkg)
   AC_GANDOLF_SETUP(pkg)
   AC_EOSPAC5_SETUP(pkg)
   AC_GRACE_SETUP(pkg)
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
dnl-------------------------------------------------------------------------dnl

