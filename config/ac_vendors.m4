dnl-------------------------------------------------------------------------dnl
dnl ac_vendors.m4
dnl macros for each vendor that is used in DRACO
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

##---------------------------------------------------------------------------##
## All vendor macros should take the following arguments:
##     pkg      - this vendor is used in the package (default)
##     test     - this vendor is used only in the package's test
##---------------------------------------------------------------------------##

dnl-------------------------------------------------------------------------dnl
dnl AC_MPI_SETUP
dnl
dnl MPI implementation (off by default)
dnl MPI is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_MPI_SETUP, [dnl

   dnl define --with-mpi
   AC_ARG_WITH(mpi,
      [  --with-mpi=[vendor,mpich] 
	                  determine MPI implementation (vendor on SGI,SUN; mpich on LINUX)])

   dnl define --with-mpi-inc and --with-mpi-lib
   AC_WITH_DIR(mpi-inc, MPI_INC, \${MPI_INC_DIR},
	       [tell where MPI includes are])
   AC_WITH_DIR(mpi-lib, MPI_LIB, \${MPI_LIB_DIR},
	       [tell where MPI libs are])

   # define MPI include path
   if test -n "${MPI_INC}" ; then
       # remember that MPI_INC has the final slash
       MPI_H="\"${MPI_INC}mpi.h\""
   elif test -z "${MPI_INC}" ; then
       MPI_H="<mpi.h>"
   fi
   
   # we define MPI_H regardless of whether a PATH is set
   AC_DEFINE_UNQUOTED(MPI_H, ${MPI_H})dnl

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

   dnl we wait to set the basic MPI libraries (if it is on) until
   dnl after checking the C4 status; these functions are performed
   dnl in ac_dracoenv.m4, section SYSTEM-SPECIFIC SETUP; we do this
   dnl here because each platform has different mpi options for
   dnl vendors and mpich

   dnl note that we used to do this in a function called AC_COMM_SET;
   dnl however, there are too many platform-dependent variables 
   dnl to continue doing this; so we do all these operations in the
   dnl platform specific section of ac_dracoenv.m4
])

dnl-------------------------------------------------------------------------dnl
dnl AC_SHMEM_SETUP
dnl
dnl SHMEM implementation (off by default)
dnl SHMEM is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_SHMEM_SETUP, [dnl
	
   dnl define --enable-shmem
   AC_ARG_ENABLE(shmem,
      [  --enable-shmem          turn on shmem (off by default)])

   dnl define --with-shmem-inc and --with-shmem-lib
   AC_WITH_DIR(shmem-inc, SHMEM_INC, \${SHMEM_INC_DIR},
	       [tell where SHMEM includes are])
   AC_WITH_DIR(shmem-lib, SHMEM_LIB, \${SHMEM_LIB_DIR},
	       [tell where SHMEM libs are])

   # define SHMEM include path
   if test -n "${SHMEM_INC}" ; then
       # remember that SHMEM_INC has the final slash
       SHMEM_H="\"${SHMEM_INC}shmem.h\""
   elif test -z "${SHMEM_INC}" ; then
       SHMEM_H="<mpp/shmem.h>"
   fi
   
   # we define SHMEM_H regardless of whether a PATH is set
   AC_DEFINE_UNQUOTED(SHMEM_H, ${SHMEM_H})dnl

   # determine if this package is needed for testing or for the
   # package
   vendor_shmem=$1
	
   # if the user sets SHMEM_INC or SHMEM_LIB and enable_shmem=no
   # then turn on enable_shmem
   if test "${enable_shmem:=no}" = no ; then
       if test -n "${SHMEM_INC}" ; then
	   enable_shmem='yes'
       elif test -n "${SHMEM_LIB}" ; then
	   enable_shmem='yes'
       fi
   fi

   dnl shmem is only available on sgis in this build system (and in
   dnl general I think)
])

dnl-------------------------------------------------------------------------dnl
dnl AC_POOMA_SETUP
dnl
dnl POOMA implementation (off by default)
dnl POOMA is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_POOMA_SETUP, [dnl
	
   dnl define --enable-pooma
   AC_ARG_ENABLE(pooma,
      [  --enable-pooma          turn on pooma (off by default)])

   dnl define --with-pooma-root
   AC_WITH_DIR(pooma-root, POOMA_ROOT, \${POOMA_ROOT_DIR},
	       [tell where root to the POOMA src and libraries are])

   dnl define --with-pooma-arch
   AC_ARG_WITH(pooma-arch,
      [  --with-pooma-arch=[POOMA_ARCH]  set pooma architecture])

   # determine if this package is needed for testing or for the
   # package
   vendor_pooma=$1
	

   if test "${enable_pooma:=no}" != no ; then
      POOMA_ROOT="${with_pooma_root}"
      AC_MSG_CHECKING(for pooma root directory)
      AC_MSG_RESULT(${with_pooma_root})	

      POOMA_ARCH="${with_pooma_arch}"
      AC_MSG_CHECKING(for POOMA_ARCH)
      AC_MSG_RESULT(${POOMA_ARCH})	
      AC_VENDORLIB_SETUP(vendor_pooma, -L${with_pooma_root}/lib/${POOMA_ARCH} -lpooma)

      CXXFLAGS="${CXXFLAGS} --restrict"
      CPPFLAGS="${CPPFLAGS} -I${POOMA_ROOT}/src"
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPRNG_SETUP
dnl
dnl SPRNG LIBRARY SETUP (on by default -lfg)
dnl SPRNG is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_SPRNG_SETUP, [dnl

   dnl define --with-sprng
   AC_ARG_WITH(sprng,
      [  --with-sprng[=lib]      determine the rng lib (lfg is default)])
	
   dnl define --with-sprng-inc and --with-sprng-lib
   AC_WITH_DIR(sprng-inc, SPRNG_INC, \${SPRNG_INC_DIR},
	       [tell where SPRNG includes are])
   AC_WITH_DIR(sprng-lib, SPRNG_LIB, \${SPRNG_LIB_DIR},
	       [tell where SPRNG libraries are])

   # define SPRNG include path
   if test -n "${SPRNG_INC}" ; then
       # remember that SPRNG_INC has the final slash
       SPRNG_H="\"${SPRNG_INC}sprng.h\""
   elif test -z "${SPRNG_INC}" ; then
       SPRNG_H="<sprng.h>"
   fi
   
   # we define SPRNG_H regardless of whether a PATH is set
   AC_DEFINE_UNQUOTED(SPRNG_H, ${SPRNG_H})dnl

   # determine if this package is needed for testing or for the 
   # package
   vendor_sprng=$1

   # define standard c interface
   AC_DEFINE(__STDC__)

   # choices are with_sprng = lfg, lcg, yes, or no

   # default (sprng is no and set to lfg by default)
   if test "${with_sprng:=lfg}" = yes ; then
       with_sprng='lfg'
   fi

   # set up the libraries
   if test "${with_sprng}" != no ; then
       if test -n "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -L${SPRNG_LIB} -l${with_sprng})
       elif test -z "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -l${with_sprng})
       fi
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_PCGLIB_SETUP
dnl
dnl PCGLIB LIBRARY SETUP (on by default)
dnl PCGLIB is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_PCGLIB_SETUP, [dnl

   dnl define --enable-pcglib
   AC_ARG_ENABLE(pcglib,        
      [  --enable-pcglib         turn on pcglib])

   dnl define --with-pcglib-lib
   AC_WITH_DIR(pcglib-lib, PCGLIB_LIB, \${PCGLIB_LIB_DIR},
	       [tell where PCGLIB libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_pcglib=$1

   # set up the libraries
   if test "${enable_pcglib:=yes}" != no ; then
       if test -z "${PCGLIB_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcglib, -lpcg_f77)
       elif test -n "${PCGLIB_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcglib, -L${PCGLIB_LIB} -lpcg_f77)
       fi
   fi

   dnl note that we add some system-specific libraries for this
   dnl vendor in AC_DRACO_ENV
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
   AC_SHMEM_SETUP(pkg)
   AC_SPRNG_SETUP(pkg)
   AC_PCGLIB_SETUP(pkg)
   dnl AC_POOMA_SETUP(pkg)
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
dnl-------------------------------------------------------------------------dnl

