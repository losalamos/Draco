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
      [  --with-mpi=[vendor,mpich,lampi] 
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
   
   # if c4=mpi and with-mpi=no explicitly then 
   # define them (mpi gets set to vendor by default)
   if test "$with_c4" = mpi ; then
       if test "$with_mpi" = no ; then
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
      [  --with-aztec=[lib]      determine the aztec lib (aztec is the default)])
 
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
      [  --with-gsl=[gsl] 
                       determine GSL lib (gsl is default)])
 
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

   # if atlas is available use it's version of cblas, 
   # otherwise use the version provided by GSL
   if test "${with_lapack}" = atlas; then
       gsl_libs='-lgsl'
   else
       gsl_libs='-lgsl -lgslcblas'
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
	   AC_VENDORLIB_SETUP(vendor_gsl, -L${GSL_LIB} ${gsl_libs})
       elif test -z "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, ${gsl_libs})
       fi

       # add GSL directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GSL_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GSL_INC}"

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
      [  --with-trilinos=[lib]    determine the trilinos implementation (aztecoo is default)])
 
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
	   AC_VENDORLIB_SETUP(vendor_trilinos, -L${TRILINOS_LIB} -l${with_trilinos} -lepetra -ltriutils)
       elif test -z "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -l${with_trilinos} -lepetra -ltriutils)
       fi

       # add TRILINOS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${TRILINOS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${TRILINOS_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SCALAPACK_SETUP
dnl
dnl SCALAPACK SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SCALAPACK_SETUP], [dnl

   dnl define --with-scalapack
   AC_ARG_WITH(scalapack,
      [  --with-scalapack=[scalapack] ])
 
   dnl define --with-scalapack-lib
   AC_WITH_DIR(scalapack-lib, SCALAPACK_LIB, \${SCALAPACK_LIB_DIR},
	       [tell where SCALAPACK libraries are])

   # set default value of scalapack includes and libs
   if test "${with_scalapack:=scalapack}" = yes ; then
       with_scalapack='scalapack'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_scalapack=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_SCALAPACK_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_scalapack}" ; then

       # no includes for scalapack

       # library path
       if test -n "${SCALAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_scalapack, -L${SCALAPACK_LIB} -lscalapack)
       elif test -z "${SCALAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_scalapack, -lscalapack)
       fi

       # add SCALAPACK directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SCALAPACK_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_BLACS_SETUP
dnl
dnl BLACS SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_BLACS_SETUP], [dnl

   dnl define --with-blacs
   AC_ARG_WITH(blacs,
      [  --with-blacs=[blacs] ])
 
   dnl define --with-blacs-lib
   AC_WITH_DIR(blacs-lib, BLACS_LIB, \${BLACS_LIB_DIR},
	       [tell where BLACS libraries are])

   # set default value of blacs includes and libs
   if test "${with_blacs:=blacs}" = yes ; then
       with_blacs='blacs'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_blacs=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_BLACS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_blacs}" ; then

       # no includes for blacs

       # library path
       if test -n "${BLACS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_blacs, -L${BLACS_LIB} -lblacsF77init -lblacsCinit -lblacs -lblacsCinit -lblacs)
       elif test -z "${BLACS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_blacs, -lblacsF77init -lblacsCinit -lblacs -lblacsCinit -lblacs)
       fi

       # add BLACS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${BLACS_LIB}"

   fi

])
dnl-------------------------------------------------------------------------dnl
dnl AC_HYPRE_SETUP
dnl
dnl HYPRE SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HYPRE_SETUP], [dnl

   dnl define --with-hypre
   AC_ARG_WITH(hypre,
      [  --with-hypre=[hypre] ])
 
   dnl define --with-hypre-inc
   AC_WITH_DIR(hypre-inc, HYPRE_INC, \${HYPRE_INC_DIR},
	       [tell where HYPRE includes are])

   dnl define --with-hypre-lib
   AC_WITH_DIR(hypre-lib, HYPRE_LIB, \${HYPRE_LIB_DIR},
	       [tell where HYPRE libraries are])

   # set default value of hypre includes and libs
   if test "${with_hypre:=hypre}" = yes ; then
       with_hypre='hypre'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_hypre=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_HYPRE_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_hypre}" ; then

       # include path
       if test -n "${HYPRE_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${HYPRE_INC}"
       fi

       # library path
       if test -n "${HYPRE_LIB}" ; then

	   AC_VENDORLIB_SETUP(vendor_hypre, -L${HYPRE_LIB} -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lHYPRE_krylov -lHYPRE_utilities)

       elif test -z "${HYPRE_LIB}" ; then

	   AC_VENDORLIB_SETUP(vendor_hypre, -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lHYPRE_krylov -lHYPRE_utilities)

       fi

       # add HYPRE directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${HYPRE_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${HYPRE_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_METIS_SETUP
dnl
dnl METIS SETUP (on by default)
dnl METIS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_METIS_SETUP], [dnl

   dnl define --with-metis
   AC_ARG_WITH(metis,
      [  --with-metis=[lib]    the metis implementation])
 
   dnl define --with-metis-inc
   AC_WITH_DIR(metis-inc, METIS_INC, \${METIS_INC_DIR},
	       [tell where METIS includes are])

   dnl define --with-metis-lib
   AC_WITH_DIR(metis-lib, METIS_LIB, \${METIS_LIB_DIR},
	       [tell where METIS libraries are])

   # set default value of metis includes and libs
   if test "${with_metis:=metis}" = yes ; then
       with_metis='metis'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_metis=$1

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_METIS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_metis}" ; then

       # include path
       if test -n "${METIS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${METIS_INC}"
       fi

       # library path
       if test -n "${METIS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_metis, -L${METIS_LIB} -l${with_metis})
       elif test -z "${METIS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_metis, -l${with_metis})
       fi

       # add METIS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${METIS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${METIS_INC}"

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
dnl AC_SPICA_SETUP
dnl
dnl SPICA LIBRARY SETUP (on by default -lSpicaCSG)
dnl SPICA is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SPICA_SETUP], [dnl

   dnl define --with-spica
   AC_ARG_WITH(spica,
      [  --with-spica[=yes]                 spica is on by default])
	
   dnl define --with-spica-inc and --with-spica-lib
   AC_WITH_DIR(spica-inc, SPICA_INC, \${SPICA_INC_DIR},
	       [tell where SPICA includes are])
   AC_WITH_DIR(spica-lib, SPICA_LIB, \${SPICA_LIB_DIR},
	       [tell where SPICA libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_spica=$1

   # define variable if spica is on
   if test "${with_spica:=yes}" != no; then
       AC_DEFINE([USE_SPICA])
   fi
])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_SPICA_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_spica}"; then

       # include path
       if test -n "${SPICA_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${SPICA_INC}"
       fi
   
       # libraries
       if test -n "${SPICA_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_spica, -L${SPICA_LIB} -lSpicaCSG)
       elif test -z "${SPICA_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_spica, -lSpicaCSG)
       fi

       # add spica directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SPICA_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SPICA_INC}"

   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_XERCES_SETUP
dnl
dnl XERCES LIBRARY SETUP
dnl xerces is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_XERCES_SETUP], [dnl

   dnl define --with-xerces
   AC_ARG_WITH(xerces,
      [  --with-xerces[=lib]      determine the XERCES xml lib (xerces-c is default)])
	
   dnl define --with-xerces-inc and --with-xerces-lib
   AC_WITH_DIR(xerces-inc, XERCES_INC, \${XERCES_INC_DIR},
	       [tell where XERCES includes are])
   AC_WITH_DIR(xerces-lib, XERCES_LIB, \${XERCES_LIB_DIR},
	       [tell where XERCES libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_xerces=$1

   # default (xerces is set to xerces-c by default)
   if test "${with_xerces:=xerces-c}" = yes ; then
       with_xerces='xerces-c'
   fi
])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_XERCES_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_xerces}"; then

       # include path
       if test -n "${XERCES_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${XERCES_INC}"
       fi
   
       # libraries
       if test -n "${XERCES_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_xerces, -L${XERCES_LIB} -l${with_xerces})
       elif test -z "${XERCES_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_xerces, -l${with_xerces})
       fi

       # add xerces directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${XERCES_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${XERCES_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_HDF5_SETUP
dnl
dnl HDF5 SETUP (on by default; 'mpi' if mpi in use, else 'serial')
dnl HDF5 is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HDF5_SETUP], [dnl

   dnl define --with-hdf5
   AC_ARG_WITH(hdf5,
      [  --with-hdf5=[serial,mpi]      determine hdf5 implementation (default:  'mpi' if mpi in use, else 'serial')])
 
   dnl define --with-hdf5-inc
   AC_WITH_DIR(hdf5-inc, HDF5_INC, \${HDF5_INC_DIR},
	       [tell where HDF5 includes are])

   dnl define --with-hdf5-lib
   AC_WITH_DIR(hdf5-lib, HDF5_LIB, \${HDF5_LIB_DIR},
	       [tell where HDF5 libraries are])

   # default (mpi if mpi is in use, else serial)
   if test "${with_hdf5:=no}" = yes ; then
       if test "${with_mpi}" != no ; then
	   with_hdf5='mpi'
       else
	   with_hdf5='serial'
       fi
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_hdf5=$1

   # define variable if hdf5 is on
   if test "${with_hdf5:=yes}" != no; then
       AC_DEFINE([USE_HDF5])
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_HDF5_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_hdf5}" ; then

       # include path
       if test -n "${HDF5_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${HDF5_INC}"
       fi

       # library path
       if test -n "${HDF5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_hdf5, -L${HDF5_LIB} -lhdf5 -lz)
       elif test -z "${HDF5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_hdf5, -lhdf5 -lz)
       fi

       # add HDF5 directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${HDF5_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${HDF5_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_UDM_SETUP
dnl
dnl UDM SETUP (on by default; 'mpi' if mpi in use, else 'serial')
dnl UDM is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_UDM_SETUP], [dnl

   dnl define --with-udm
   AC_ARG_WITH(udm,
      [  --with-udm=[serial,mpi]      determine udm implementation (default:  'mpi' if mpi in use, else 'serial')])
 
   dnl define --with-udm-inc
   AC_WITH_DIR(udm-inc, UDM_INC, \${UDM_INC_DIR},
	       [tell where UDM includes are])

   dnl define --with-udm-lib
   AC_WITH_DIR(udm-lib, UDM_LIB, \${UDM_LIB_DIR},
	       [tell where UDM libraries are])

   # default (mpi if mpi is in use, else serial)
   if test "${with_udm:=no}" = yes ; then
       if test "${with_mpi}" != no ; then
	   with_udm='mpi'
       else
	   with_udm='serial'
       fi
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_udm=$1

   # define variable if udm is on
   if test "${with_udm:=no}" != no; then
       AC_DEFINE([USE_UDM])
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_UDM_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_udm}" ; then

       # include path
       if test -n "${UDM_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${UDM_INC}"
           # set extra #define if using udm in parallel
           if test "${with_udm}" = mpi ; then
               AC_DEFINE(UDM_HAVE_PARALLEL)
           fi
       fi

       # library path
       if test -n "${UDM_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_udm, -L${UDM_LIB} -ludm)
       elif test -z "${UDM_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_udm, -ludm)
       fi

       # add UDM directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${UDM_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${UDM_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DLOPEN_SETUP
dnl
dnl This is an optional vendor.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DLOPEN_SETUP], [dnl

   dnl define --enable-dlopen
   AC_ARG_ENABLE(dlopen,
      [  --enable-dlopen          Enable dlopen (default: on if --enable-shared, off otherwise)])

   # determine if this package is needed for testing or for the
   # package.
   vendor_dlopen=$1 

   # set default value for enable_dlopen, which is the value of enable_shared.
   if test "${enable_shared}" = yes ; then
       if test "${enable_dlopen:=yes}" != no ; then 
	   enable_dlopen=yes
       fi
   else
       if test "${enable_dlopen:=no}" != no ; then 
	   enable_dlopen=yes
       fi
   fi

   # turn off dlopen if not using shared libraries.
   if test "${enable_shared}" != yes ; then
       if test "${enable_dlopen}" = yes ; then
	   AC_MSG_WARN("Must specify --enable-shared when using --enable-dlopen.")
           AC_MSG_WARN("   dlopen disabled.")
       fi
       enable_dlopen=no
   fi

   if test "${enable_dlopen}" = yes ; then
       AC_DEFINE(USE_DLOPEN)
   fi
]) 

##---------------------------------------------------------------------------##

AC_DEFUN([AC_DLOPEN_FINALIZE], [dnl
   # Libraries are platform-specific; done in ac_platforms.
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
   AC_GSL_FINALIZE

   AC_AZTEC_FINALIZE
   AC_PCG_FINALIZE
   AC_HYPRE_FINALIZE
   AC_SCALAPACK_FINALIZE
   AC_BLACS_FINALIZE
   AC_LAPACK_FINALIZE
   AC_EOSPAC5_FINALIZE
   AC_GANDOLF_FINALIZE
   AC_SPRNG_FINALIZE
   AC_GRACE_FINALIZE
   AC_METIS_FINALIZE
   AC_SPICA_FINALIZE
   AC_XERCES_FINALIZE

   AC_UDM_FINALIZE
   AC_HDF5_FINALIZE

   AC_MPI_FINALIZE
   AC_DLOPEN_FINALIZE

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
   AC_TRILINOS_SETUP(pkg)
   AC_METIS_SETUP(pkg)
   AC_LAPACK_SETUP(pkg)
   AC_GANDOLF_SETUP(pkg)
   AC_EOSPAC5_SETUP(pkg)
   AC_GRACE_SETUP(pkg)
   AC_SPICA_SETUP(pkg)
   AC_XERCES_SETUP(pkg)
   AC_HDF5_SETUP(pkg)
   AC_UDM_SETUP(pkg)
   AC_DLOPEN_SETUP(pkg)
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
dnl-------------------------------------------------------------------------dnl

