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
dnl MPI implementation (on by default)
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

   # make sure that the host is defined
   AC_REQUIRE([AC_CANONICAL_HOST])

   # platform defaults:
   if test ${with_mpi:-no} = no || test ${with_mpi:-no} = yes; then
   case $host in
   mips-sgi-irix6.*) 
      mpi_default=vendor ;;
   alpha*-dec-osf*)
      mpi_default=vendor ;;
   *ibm-aix*)
      mpi_default=vendor ;;
   *-linux-gnu)
      mpi_default=openmpi ;;
   *)
      mpi_default=vendor ;;
   esac
   else
     mpi_default=${with_mpi}
   fi

   dnl echo "Before: with_mpi = ${with_mpi}"

   AC_SETUP_VENDOR( [mpi], [yes], [${mpi_default}], 
                    [@S|@{MPI_INC_DIR}], [@S|@{MPI_LIB_DIR}], 
                    [vendor|mpich|lampi|openmpi] )

   dnl echo "After: with_mpi = ${with_mpi}"
   dnl echo "with_c4  = ${with_c4}"

   if test "${with_mpi}" = "no" ; then
      # if with_mpi is no, but with_c4=mpi, then abort.
      if test "${with_c4}" = "mpi"; then
        { echo "configure: error: --with-c4=mpi requires with_mpi to \
have a value other than 'no'.  Try setting --with-c4=scalar or \
setting --with-mpi=[vendor|mpich|lampi|openmpi].  Also examine the \
values of @S|@MPI_INC_DIR and @S|@MPI_LIB_DIR." 1>&2; \
          exit 1; }
      fi
   else
      # determine if this package is needed for testing or for the
      # package
      vendor_mpi=$1 
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_MPI_FINALIZE], [dnl

   dnl echo "mpi_finalize: with_mpi=${with_mpi}"
   dnl echo "              with_c4 =${with_c4}"

   save_with_mpi=${with_mpi}
   case ${with_mpi} in
   mpich) 
     dnl default works correctly.
     ;;
   lampi | LAMPI | openmpi | OPENMPI)
      with_mpi=mpi dnl We want '-lmpi' on the link line
      ;;
   esac

   dnl mpi_libs is set in ac_vendors.m4
   AC_FINALIZE_VENDOR([mpi],[${mpi_libs}]) 

   dnl reset with_mpi value.
   with_mpi=${save_with_mpi}

])

dnl-------------------------------------------------------------------------dnl
dnl AC_NDI_SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_NDI_SETUP], [dnl

   dnl define --with-ndi
   AC_ARG_WITH([ndi],
     [AS_HELP_STRING([--with-ndi@<:@=ndi@:>@],
       [determine NDI lib (ndi is default)])])

   dnl define --with-ndi-inc
   AC_WITH_DIR(ndi-inc, NDI_INC, @S|@{NDI_INC_DIR},
               [tell where NDI includes are])

   dnl define --with-ndi-lib
   AC_WITH_DIR(ndi-lib, NDI_LIB, @S|@{NDI_LIB_DIR},
               [tell where NDI libraries are])

   # set default value of ndi includes and libs
   if test "${with_ndi:=ndi}" = yes ; then
       with_ndi='ndi'
   fi

   ndi_libs='-lndi'

   # determine if this package is needed for testing or for the package
   vendor_ndi=$1
])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_NDI_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_ndi}"; then

       # include path
       if test -n "${NDI_INC}"; then
           # add to include path
           VENDOR_INC="${VENDOR_INC} -I${NDI_INC}"
       fi

       # library path
       if test -n "${NDI_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_ndi, -L${NDI_LIB} ${ndi_libs})
       elif test -z "${NDI_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_ndi, ${ndi_libs})
       fi

       # add NDI directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${NDI_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${NDI_INC}"
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
   AC_ARG_WITH([aztec],
     [AS_HELP_STRING([--with-aztec@<:@=aztec@:>@],
       [determine the aztec lib (aztec is the default)])])
 
   dnl define --with-aztec-inc
   AC_WITH_DIR(aztec-inc, AZTEC_INC, @S|@{AZTEC_INC_DIR},
               [tell where AZTEC includes are])

   dnl define --with-aztec-lib
   AC_WITH_DIR(aztec-lib, AZTEC_LIB, @S|@{AZTEC_LIB_DIR},
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

AC_DEFUN([AC_GSL_SETUP], [dnl

   AC_SETUP_VENDOR( gsl, yes, gsl,
                   @S|@{GSL_INC_DIR}, @S|@{GSL_LIB_DIR},
                   [gsl] )

   # Set up gsl only if --with-gsl or --with-gsl-lib is
   # explicitly set if $with_gsl is "yes" or if it has a value other
   # than "no" the setup gsl.  $with_gsl will be "no" if
   # --without-gsl is specified on the configure line.
   if ! test "${with_gsl}" = "no" ; then

     # if atlas is available use it's version of cblas, 
     # otherwise use the version provided by GSL
     if ! test "${with_lapack}" = atlas; then
       gsl_extra_libs='-lgslcblas'
     fi

     # determine if this package is needed for testing or for the 
     # package
     vendor_gsl=$1

   fi
])

##---------------------------------------------------------------------------##

dnl AC_DEFUN([AC_GSL_FINALIZE], [dnl
dnl    AC_FINALIZE_VENDOR([gsl],[${gsl_extra_libs}])
dnl ])

dnl-------------------------------------------------------------------------dnl
dnl AC_SUPERLUDIST_SETUP
dnl
dnl SUPERLUDIST SETUP (on by default)
dnl SUPERLUDIST is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SUPERLUDIST_SETUP], [dnl

   dnl define --with-superludist
   AC_ARG_WITH([superludist],
     [AS_HELP_STRING([--with-superludist@<:@=superludist@:>@],
       [determine SUPERLUDIST lib (superludist is default)])])

   dnl define --with-superludist-inc
   AC_WITH_DIR(superludist-inc, SUPERLUDIST_INC, @S|@{SUPERLUDIST_INC_DIR},
               [tell where SUPERLUDIST includes are])

   dnl define --with-superludist-lib
   AC_WITH_DIR(superludist-lib, SUPERLUDIST_LIB, @S|@{SUPERLUDIST_LIB_DIR},
               [tell where SUPERLUDIST libraries are])

   # set default value of superludist includes and libs
   if test "${with_superludist:=superludist}" = yes ; then
      with_superludist='superludist'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_superludist=$1
])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_SUPERLUDIST_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_superludist}"; then

       # include path
       if test -n "${SUPERLUDIST_INC}"; then 
           # add to include path
           VENDOR_INC="${VENDOR_INC} -I${SUPERLUDIST_INC}"
       fi

       # library path
       # if this is a scalar build, use SUPERLU instead.
       if test "${with_c4}" = "scalar" ; then
         if test -n "${SUPERLUDIST_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_superludist, -L${SUPERLUDIST_LIB} -lsuperlu)
         elif test -z "${SUPERLUDIST_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_superludist, -lsuperlu)
         fi
       else
         if test -n "${SUPERLUDIST_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_superludist, -L${SUPERLUDIST_LIB} -lsuperludist)
         elif test -z "${SUPERLUDIST_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_superludist, -lsuperludist)
         fi
       fi

       # add SUPERLUDIST directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SUPERLUDIST_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SUPERLUDIST_INC}"

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
   AC_ARG_WITH([trilinos],
     [AS_HELP_STRING([--with-trilinos@<:@=aztecoo@:>@],
       [determine the trilinos implementation (aztecoo is default)])])
 
   dnl define --with-trilinos-inc
   AC_WITH_DIR(trilinos-inc, TRILINOS_INC, @S|@{TRILINOS_INC_DIR},
               [tell where TRILINOS includes are])

   dnl define --with-trilinos-lib
   AC_WITH_DIR(trilinos-lib, TRILINOS_LIB, @S|@{TRILINOS_LIB_DIR},
               [tell where TRILINOS libraries are])

   # set default value of trilinos includes and libs
   if test "${with_trilinos:=yes}" = yes ; then
       with_trilinos='-llocathyra -llocaepetra -lloca -lnoxthyra -lnoxepetra -lnox -lModeLaplace -lanasaziepetra -lanasazi -lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml -lstratimikosifpack -lbelosepetra -lbelos -lml -lifpack -lamesos -lgaleri -laztecoo -lthyraepetraext -lthyraepetra -lthyra -lepetraext -ltriutils -lzoltan -lepetra -lrtop -lteuchos'  
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
           AC_VENDORLIB_SETUP(vendor_trilinos, -L${TRILINOS_LIB} ${with_trilinos})
       elif test -z "${TRILINOS_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_trilinos, ${with_trilinos})
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
   AC_ARG_WITH([scalapack],
     [AS_HELP_STRING([--with-scalapack@<:@=scalapack@:>@],
       [determine the scalapack library name (default is scalapack)])])
 
   dnl define --with-scalapack-lib
   AC_WITH_DIR(scalapack-lib, SCALAPACK_LIB, @S|@{SCALAPACK_LIB_DIR},
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
   AC_ARG_WITH([blacs],[  --with-blacs=[blacs] ])
 
   dnl define --with-blacs-lib
   AC_WITH_DIR(blacs-lib, BLACS_LIB, @S|@{BLACS_LIB_DIR},
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
   AC_ARG_WITH([hypre],[  --with-hypre=[hypre] ])
 
   dnl define --with-hypre-inc
   AC_WITH_DIR(hypre-inc, HYPRE_INC, @S|@{HYPRE_INC_DIR},
               [tell where HYPRE includes are])

   dnl define --with-hypre-lib
   AC_WITH_DIR(hypre-lib, HYPRE_LIB, @S|@{HYPRE_LIB_DIR},
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

           AC_VENDORLIB_SETUP(vendor_hypre, -L${HYPRE_LIB} -lHYPRE)

       elif test -z "${HYPRE_LIB}" ; then

           AC_VENDORLIB_SETUP(vendor_hypre, -lHYPRE)

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
   AC_ARG_WITH([metis],[  --with-metis=[metis]      the metis implementation])
 
   dnl define --with-metis-inc
   AC_WITH_DIR(metis-inc, METIS_INC, @S|@{METIS_INC_DIR},
               [tell where METIS includes are])

   dnl define --with-metis-lib
   AC_WITH_DIR(metis-lib, METIS_LIB, @S|@{METIS_LIB_DIR},
               [tell where METIS libraries are])

   # set default value of metis includes and libs
   if test "${with_metis:=metis}" = yes ; then
       with_metis='metis'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_metis=$1

])

dnl-------------------------------------------------------------------------dnl
dnl AC_PARMETIS_SETUP
dnl
dnl PARMETIS SETUP (on by default)
dnl PARMETIS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_PARMETIS_SETUP], [dnl

   dnl define --with-parmetis
   AC_ARG_WITH([parmetis],
     [AS_HELP_STRING([--with-parmetis@<:@=parmetis@:>@],
       [the parmetis implementation])])
 
   dnl define --with-parmetis-inc
   AC_WITH_DIR(parmetis-inc, PARMETIS_INC, @S|@{PARMETIS_INC_DIR},
               [tell where PARMETIS includes are])

   dnl define --with-parmetis-lib
   AC_WITH_DIR(parmetis-lib, PARMETIS_LIB, @S|@{PARMETIS_LIB_DIR},
               [tell where PARMETIS libraries are])

   # set default value of parmetis includes and libs
   if test "${with_parmetis:=parmetis}" = yes ; then
       with_parmetis='parmetis'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_parmetis=$1

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

##---------------------------------------------------------------------------##

AC_DEFUN([AC_PARMETIS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_parmetis}" ; then

       # include path
       if test -n "${PARMETIS_INC}"; then 
           # add to include path
           VENDOR_INC="${VENDOR_INC} -I${PARMETIS_INC}"
       fi

       # library path
       if test -n "${PARMETIS_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_parmetis, -L${PARMETIS_LIB} -l${with_parmetis})
       elif test -z "${PARMETIS_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_parmetis, -l${with_parmetis})
       fi

       # add PARMETIS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${PARMETIS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${PARMETIS_INC}"

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

   AC_SETUP_VENDOR( [pcg], [no], [pcg],
                   [], [@S|@{PCG_LIB_DIR}], [pcg] )

   # Set up pcg only if --with-pcg or --with-pcg-lib is explicitly set
   if ! test  "${with_pcg:-no}" = no ; then

      AC_DEFINE(USE_PCGLIB)

      # determine if this package is needed for testing or for the 
      # package
      vendor_pcg=$1
   fi

])

##---------------------------------------------------------------------------##

dnl AC_DEFUN([AC_PCG_FINALIZE], [dnl
dnl    AC_FINALIZE_VENDOR([pcg],[${pcg_extra_libs}])
dnl ])

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

   AC_SETUP_VENDOR( [gandolf], [yes], [gandolf],
                   [], [@S|@{GANDOLF_LIB_DIR}],
                   [gandolf] )

   # gandolf is set to libgandolf by default
   dnl echo "gandolf_setup: with_gandolf = ${with_gandolf}"

   if ! test "${with_gandolf:=gandolf}" = "no" ; then
       # determine if this package is needed for testing or for the 
       # package
       vendor_gandolf=$1
   else
       ## TK added May 7 07: 
       vendor_gandolf=''                   ## skip finalize
       with_gandolf=''                     ## stub this out
       AC_MSG_RESULT("NOT USING GANDOLF")  ## alert/remind user
       AC_DEFINE(rtt_cdi_gandolf_stub)     ## used in cdi_gandolf/config.h
   fi
])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_GANDOLF_FINALIZE], [dnl

   # set up the libraries
   # This string should be blank, 'test' or 'pkg'
   # This is set except when $with_gandolf = 'no'
   if test -n "${vendor_gandolf}" ; then

     # [2010-04-28 KT] Special case when on Linux and libgfortran is
     # available, but is incompatible with libgandolf.a.  This is the
     # case on ccscs8 when using gcc-4.3.4 while libgandolf.a expects
     # gcc-4.1.2 (the symbol 'gfortran_copy_string' is missing from
     # the newer libgfortran.so)
     case $host in
     *-linux-gnu)

       # The code in ac_platforms should have already identified $F90
dnl    if test "${F90}" = gfortran; then
dnl          # Examine libgfortran to determine if it has the features
dnl          # needed by libgandolf.a
dnl          AC_CHECK_LIB( [gfortran], [gfortran_copy_string],
dnl            [], [gandolf_gfortran_special_lib=yes] )

         # This library is crazy.  Some versions need -lg2c and others
         # need -lgfortran.  Let's try to figure it out.
         
         AC_PATH_PROG( NM_BIN, nm, null )
         AC_MSG_CHECKING([for extra libraries to support vendor gandolf])
         if test -x ${NM_BIN}; then
           if test -n "${with_gandolf_lib}"; then
             # undefined symbols found in libgandolf.a
             libgandolf_need_gfortran=`$NM_BIN -a ${with_gandolf_lib}/libgandolf.a | grep " U _gfortran_compare_string"`
             libgandolf_need_g2c=`$NM_BIN -a ${with_gandolf_lib}/libgandolf.a | grep " U s_copy"`
             if test -n "${libgandolf_need_gfortran}"; then
               gandolf_gfortran_special_lib='-lgfortran'
             fi
             if test -n "${libgandolf_need_g2c}"; then
               gandolf_gfortran_special_lib='${gandolf_gfortran_special_lib} -lg2c'
             fi
           fi
         fi     

         AC_MSG_RESULT([${gandolf_gfortran_special_lib}])
       ;;
     esac
   fi

   AC_FINALIZE_VENDOR([gandolf],[${gandolf_gfortran_special_lib}])
])

dnl-------------------------------------------------------------------------dnl
dnl AC_EOSPAC5_SETUP
dnl
dnl EOSPAC5 LIBRARY SETUP (on by default)
dnl EOSPAC5 is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_EOSPAC5_SETUP], [dnl

   dnl define --with-eospac
   AC_ARG_WITH([eospac],
     [AS_HELP_STRING([--with-eospac@<:@=eospac@:>@],
       [determine the eospac lib name (eospac is default)])])

   dnl define --with-eospac-lib
   AC_WITH_DIR(eospac-lib, EOSPAC5_LIB, @S|@{EOSPAC5_LIB_DIR},
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

   # By default use vendor supplied LAPACK.  On Linux, if
   # LAPACK_LIB_DIR is set to an atlas path, then use atlas by
   # default. 
   lapack_default="vendor"
   if test -n "`echo $LAPACK_LIB_DIR | grep atlas`"; then
      lapack_default="atlas"
   fi
   if test ${with_lapack} = atlas; then
      lapack_default=$with_lapack
   fi

   AC_SETUP_VENDOR( [lapack], [yes], [${lapack_default}],
                   [], [@S|@{LAPACK_LIB_DIR}], [vendor|atlas] )

   if ! test "${with_lapack}" = no ; then

      # define the atlas libraries (these are system independent)
      if test "${with_lapack}" = atlas; then
         lapack_extra_libs='-llapack -lf77blas -lcblas -latlas'
      fi

      # If lapack is built with pgf90, then we may need more libraries.
#      if test "${with_cxx}" = pgi; then
#         lapack_extra_libs="${lapack_extra_libs} -lpgf90"
#      fi

      # determine if this package is needed for testing or for the
      # package
      vendor_lapack=$1

   fi

])

##---------------------------------------------------------------------------##

dnl AC_DEFUN([AC_LAPACK_FINALIZE], [dnl
dnl    AC_FINALIZE_VENDOR([lapack],[${lapack_extra_libs}])
dnl ])

dnl-------------------------------------------------------------------------dnl
dnl AC_GRACE_SETUP
dnl
dnl GRACE SETUP (on by default)
dnl GRACE is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GRACE_SETUP], [dnl

   AC_SETUP_VENDOR( grace, yes, grace_np, 
                   @S|@{GRACE_INC_DIR}, @S|@{GRACE_LIB_DIR}, [grace_np] )

   if ! test "${with_grace}" = "no" ; then

      # define GRACE header file
      GRACE_H="<${with_grace}.h>"
      AC_DEFINE_UNQUOTED(GRACE_H, ${GRACE_H})dnl

      # determine if this package is needed for testing or for the 
      # package
      vendor_grace=$1
   fi

])

##---------------------------------------------------------------------------##

dnl AC_DEFUN([AC_GRACE_FINALIZE], [dnl
dnl    AC_FINALIZE_VENDOR( grace )
dnl ])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPICA_SETUP
dnl
dnl SPICA LIBRARY SETUP (on by default -lSpicaCSG)
dnl SPICA is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SPICA_SETUP], [dnl

   dnl define --with-spica
   AC_ARG_WITH([spica],
     [AS_HELP_STRING([--with-spica],[spica is on by default])])
        
   dnl define --with-spica-inc and --with-spica-lib
   AC_WITH_DIR(spica-inc, SPICA_INC, @S|@{SPICA_INC_DIR},
               [tell where SPICA includes are])
   AC_WITH_DIR(spica-lib, SPICA_LIB, @S|@{SPICA_LIB_DIR},
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
   AC_ARG_WITH([xerces],
     [AS_HELP_STRING([--with-xerces@<:@=xerces-c@:>@],
       [determine the XERCES xml lib (xerces-c is default)])])
        
   dnl define --with-xerces-inc and --with-xerces-lib
   AC_WITH_DIR(xerces-inc, XERCES_INC, @S|@{XERCES_INC_DIR},
               [tell where XERCES includes are])
   AC_WITH_DIR(xerces-lib, XERCES_LIB, @S|@{XERCES_LIB_DIR},
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
   AC_ARG_WITH([hdf5],
     [AS_HELP_STRING([--with-hdf5@<:@=serial|mpi@:>@],
       [determine hdf5 implementation (default is mpi if mpi in use, otherwise serial)])])
 
   dnl define --with-hdf5-inc
   AC_WITH_DIR(hdf5-inc, HDF5_INC, @S|@{HDF5_INC_DIR},
               [tell where HDF5 includes are])

   dnl define --with-hdf5-lib
   AC_WITH_DIR(hdf5-lib, HDF5_LIB, @S|@{HDF5_LIB_DIR},
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
   AC_ARG_WITH([udm],
     [AS_HELP_STRING([--with-udm@<:@=serial|mpi@:>@],
       [determine udm implementation (default is mpi if mpi in use, else serial)])])
 
   dnl define --with-udm-inc
   AC_WITH_DIR(udm-inc, UDM_INC, @S|@{UDM_INC_DIR},
               [tell where UDM includes are])

   dnl define --with-udm-lib
   AC_WITH_DIR(udm-lib, UDM_LIB, @S|@{UDM_LIB_DIR},
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
dnl AC_SILO_SETUP
dnl
dnl SILO SETUP
dnl SILO is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SILO_SETUP], [dnl

   dnl define --with-silo
   AC_ARG_WITH([silo],
     [AS_HELP_STRING([--with-silo],
       [use silo (default is on)])])

   dnl define --with-silo-inc
   AC_WITH_DIR(silo-inc, SILO_INC, @S|@{SILO_INC_DIR},
               [tell where SILO includes are])

   dnl define --with-silo-lib
   AC_WITH_DIR(silo-lib, SILO_LIB, @S|@{SILO_LIB_DIR},
               [tell where SILO libraries are])

   dnl if either silo-inc or silo-lib defined, then set with_silo
   dnl thus, don't need --with-silo when using --with-silo-inc or
   dnl --with-silo-lib
   if test -n "${SILO_INC}" ; then
    with_silo="yes"
   fi
   if test -n "${SILO_LIB}" ; then
    with_silo="yes"
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_silo=$1

   # define variable if silo is on
   if test "${with_silo:=no}" != no; then
       AC_DEFINE([USE_SILO])
   fi

])

##---------------------------------------------------------------------------##

AC_DEFUN([AC_SILO_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_silo}" ; then

       # include path
       if test -n "${SILO_INC}"; then
           # add to include path
           VENDOR_INC="${VENDOR_INC} -I${SILO_INC}"
       fi

       # library path
       if test -n "${SILO_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_silo, -L${SILO_LIB} -lsiloh5)
       elif test -z "${SILO_LIB}" ; then
           AC_VENDORLIB_SETUP(vendor_silo, -lsiloh5)
       fi

       # add SILO directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SILO_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SILO_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DLOPEN_SETUP
dnl
dnl This is an optional vendor.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DLOPEN_SETUP], [dnl

   dnl define --enable-dlopen
   AC_ARG_ENABLE([dlopen],
     [AS_HELP_STRING([--enable-dlopen],
       [Enable dlopen (default: on if --enable-shared, off otherwise)])])

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

   AC_TRILINOS_FINALIZE dnl  Depends on: LAPACK, MPI
   dnl AC_GSL_FINALIZE           dnl Depends on: LAPACK
   AC_FINALIZE_VENDOR([gsl],[${gsl_extra_libs}])dnl Depends on: LAPACK

   AC_AZTEC_FINALIZE

   AC_SUPERLUDIST_FINALIZE   dnl Depends on: PARMETIS
   AC_PARMETIS_FINALIZE      dnl Depends on: METIS
   AC_METIS_FINALIZE

   dnl AC_PCG_FINALIZE        dnl Depends on: LAPACK
   AC_FINALIZE_VENDOR([pcg])  dnl Depends on: LAPACK

   AC_HYPRE_FINALIZE
   AC_SCALAPACK_FINALIZE     dnl Depends on: BLACS, MPI
   AC_BLACS_FINALIZE         dnl Depends on: MPI

   dnl AC_LAPACK_FINALIZE
   AC_FINALIZE_VENDOR([lapack],[${lapack_extra_libs}]) dnl Optionally depends on: ATLAS
   AC_EOSPAC5_FINALIZE
   AC_GANDOLF_FINALIZE
   dnl AC_GRACE_FINALIZE
   AC_FINALIZE_VENDOR([grace])
   AC_SPICA_FINALIZE
   AC_XERCES_FINALIZE

   AC_NDI_FINALIZE

   AC_UDM_FINALIZE
   AC_SILO_FINALIZE
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

AC_DEFUN([AC_ALL_VENDORS_SETUP], [dnl

   dnl include all macros for easy use in top-level configure.in's
   AC_MPI_SETUP(pkg)
   AC_PCG_SETUP(pkg)
   AC_AZTEC_SETUP(pkg)
   AC_GSL_SETUP(pkg)
   AC_SUPERLUDIST_SETUP(pkg)
   AC_TRILINOS_SETUP(pkg)
   AC_PARMETIS_SETUP(pkg)
   AC_METIS_SETUP(pkg)
   AC_LAPACK_SETUP(pkg)
   AC_GANDOLF_SETUP(pkg)
   AC_EOSPAC5_SETUP(pkg)
   AC_NDI_SETUP(pkg)
   AC_GRACE_SETUP(pkg)
   AC_SPICA_SETUP(pkg)
   AC_XERCES_SETUP(pkg)
   AC_HDF5_SETUP(pkg)
   AC_UDM_SETUP(pkg)
   AC_SILO_SETUP(pkg)
   AC_DLOPEN_SETUP(pkg)
  
])

dnl-------------------------------------------------------------------------dnl
dnl AC_VENDORS_REPORT
dnl
dnl DRACO Report of vendors found
dnl-------------------------------------------------------------------------dnl
AC_DEFUN([AC_VENDOR_REPORT], [dnl

   echo " "
   echo "Configuration Report:"
   echo " "
   echo "   Prefix                 : ${prefix}"
   echo "   Debug symbols          : ${enable_debug}"
   echo "   Optimization level     : ${with_opt}"
   echo "   Design-by-Contract     : ${with_dbc}"
dnl this is not defined in the top level configure.ac
dnl echo "   Compiler vendor        : ${with_cxx}" 
   echo "   Create shared libraries: ${enable_shared}"

   echo " "
   echo "Vendors found:"
   echo " "
   if test ${with_mpi:-no} != no; then
      echo "   MPI                YES - ${with_mpi}"
   else
      echo "   MPI                NO  - all code will be built with-c4=scalar."
   fi
   if test ${with_pcg:-no} != no; then
      echo "   PCG                YES"
   else
      echo "   PCG                NO  - pcgWrap will be omitted."
   fi
   if test ${with_gsl:-no} != no; then
      echo "   GSL                YES"
   else
      echo "   GSL                NO  - special_functions and quadrature will be omitted."
   fi
   if test ${with_lapack:-no} != no; then
      echo "   LAPACK             YES"
   else
      echo "   LAPACK             NO  - lapack_wrap will be omitted."
   fi
   if test ${with_gandolf:-no} != no; then
      echo "   GANDOLF            YES"
   else
      echo "   GANDOLF            NO  - cdi_gandolf will be omitted."
   fi
   if test ${with_grace:-no} != no; then
      echo "   GRACE              YES"
   else
      echo "   GRACE              NO  - plot2D will be omitted."
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
dnl-------------------------------------------------------------------------dnl

