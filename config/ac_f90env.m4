dnl ========================================================================
dnl 
dnl 	Author:	Mark G. Gray
dnl 		Los Alamos National Laboratory
dnl 	Date:	Wed Apr 19 16:39:19 MDT 2000
dnl 
dnl 	Copyright (c) 2000 U. S. Department of Energy. All rights reserved.
dnl 
dnl ========================================================================

dnl NAME

dnl	AC_WITH_F90, AC_F90_ENV

dnl SYNOPSIS/USAGE

dnl     AC_WITH_F90
dnl     AC_F90_ENV

dnl DESCRIPTION

dnl     AC_WITH_F90 sets the variable with_f90 to yes if it is not already 
dnl     set.

dnl     AC_F90_ENV set environment variables F90, F90FLAGS, F90EXT, 
dnl     F90FREE, F90FIXED, and MODFLAG for the compiler requested by 
dnl     with_f90.  If no specific compiler is requested, guess a compiler 
dnl     based on the target
dnl ========================================================================

##---------------------------------------------------------------------------##
## $Id$ 
##---------------------------------------------------------------------------##

dnl ### Ensure with_f90 set
AC_DEFUN([AC_WITH_F90], [dnl
   : ${with_f90:=yes}
    
   dnl turn off C++ compiler
dnl   with_cxx='no'

   dnl defines --with-f90

   AC_ARG_WITH([f90],
     [AS_HELP_STRING([--with-f90@<:@string@:>@],
      [choose a F90 compiler. Currently accepted values are
       @<:@ gfortran | XL | Fujitsu | Lahey | Portland | WorkShop |
       Cray | MIPS | Compaq | HP | Intel | NAG | Absoft @:>@ ])])
])

dnl KT (2010-04-27): We should probably inspect the value of ${FC} if
dnl --with-f90=yes is provided.

dnl
dnl CHOOSE A F90 COMPILER
dnl

AC_DEFUN([AC_F90_ENV], [dnl
   AC_REQUIRE([AC_CANONICAL_HOST])

   case "${with_f90:=yes}" in
   XL)
       AC_COMPILER_XL_F90
   ;;
   Fujitsu)
       AC_COMPILER_FUJITSU_F90
   ;;
   gfortran)
       AC_COMPILER_GFORTRAN_F90
   ;;
   Lahey)
       AC_COMPILER_LAHEY_F90
   ;;
   Portland)
       AC_COMPILER_PORTLAND_F90
   ;;
   WorkShop)
       AC_COMPILER_WORKSHOP_F90
   ;;
   Cray)
      AC_COMPILER_CRAY_F90
   ;;
   MIPS)
       AC_COMPILER_MIPS_F90
   ;;
   Compaq)
       AC_COMPILER_COMPAQ_F90
   ;;
   HP)
       AC_COMPILER_HP_F90
   ;;
   Intel)
       AC_COMPILER_INTEL_F90
   ;;
   NAG)
       AC_COMPILER_NAG_F90
   ;;
   Absoft)
       AC_COMPILER_ABSOFT_F90
   ;;
   yes)				# guess compiler from target platform
       case "${host}" in   
       rs6000-ibm-aix*)
           AC_COMPILER_XL_F90
       ;;
       powerpc-ibm-aix*)
           AC_COMPILER_XL_F90
       ;;
       sparc-sun-solaris2.*)
           AC_COMPILER_WORKSHOP_F90
       ;;
       i?86-pc-linux*)
           AC_COMPILER_LAHEY_F90
       ;;
       ymp-cray-unicos*)
          AC_COMPILER_CRAY_F90
       ;;
       mips-sgi-irix*)
          AC_COMPILER_MIPS_F90
       ;;
       i??86-pc-cygwin*)
          AC_COMPILER_COMPAQ_F90
       ;;
       alpha*)
          AC_COMPILER_COMPAQ_F90
       ;;
       *hp-hpux*)
          AC_COMPILER_HP_F90
       ;;
       *)
          AC_MSG_ERROR([Cannot guess F90 compiler, set --with-f90])
       ;;
       esac
   ;;
   no)
   ;;
   *)
       AC_MSG_ERROR([Unrecognized F90 compiler, use --help])
   ;;
   esac

   # Always use the CXX linker if the project has C++.  If this is a
   # F90 only project, then use the F90 linker.
   if test ${with_cxx:-yes} = no; then
      AR=${F90AR}
      LD=${F90LD}
   fi

   AC_SUBST(F90FREE)
   AC_SUBST(F90FIXED)
   AC_SUBST(F90FLAGS)
   AC_SUBST(MODFLAG)
])

dnl-------------------------------------------------------------------------dnl
dnl IBM XLF95 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_XL_F90], [dnl

   # Check for working XL F90 compiler

  if test "${with_upslib:=no}" != "no"
  then
     AC_CHECK_PROG(F90, mpxlf95, mpxlf95, none)
     if test "${F90}" != mpxlf95
     then
         AC_MSG_ERROR([not found])
     fi
  else
     AC_CHECK_PROG(F90, xlf95, xlf95, none)
     if test "${F90}" != xlf95
     then
         AC_MSG_ERROR([not found])
     fi
  fi
  
   # FREE, FIXED AND MODULE FLAGS

   F90FREE='-qfree=f90'
   F90FIXED='-qfixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="-qsuffix=f=f90 -qmaxmem=-1 -qextchk -qarch=pwr2 -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp -qhalt=s ${F90FREE}"
       F90FLAGS="-qsuffix=f=f90 -qmaxmem=-1 -qextchk -qarch=auto -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp -qnosave -qlanglvl=95pure -qzerosize ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	   trapflags="-qinitauto=FF"
	   trapflags="${trapflags} -qflttrap=overflow:underflow:zerodivide:invalid:enable"
	   trapflags="${trapflags} -qsigtrap"
	   F90FLAGS="-g -d -C ${trapflags} -bloadmap:loadmap.dat ${F90FLAGS}"
       else
	 # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	   F90FLAGS="-O3 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_XL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl FUJITSU F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_FUJITSU_F90], [dnl

   # Check for working Fujitsu F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "Fujitsu"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Free'
   F90FIXED='-Fixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="-X9 -Am ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -Haesu ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_FUJITSU_F90
])

dnl-------------------------------------------------------------------------dnl
dnl LAHEY F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_LAHEY_F90], [dnl

   AC_CHECK_PROG(F90, lf95, lf95, none)

   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='--nfix'
   F90FIXED='--fix'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="--f95 ${F90FREE}"
       F90FLAGS="--staticlink --f95 --in --info --swm 2004,2006,2008,8202,8203,8204,8205,8206,8209,8220 ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	  # F90FLAGS="-g --chk --trace ${F90FLAGS}"
	    F90FLAGS="-g --ap --chk --pca --private --trap --wo ${F90FLAGS}"
       else
	  # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	    F90FLAGS="-O --ntrace ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_LAHEY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl GFORTRAN F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_GFORTRAN_F90], [dnl

   AC_MSG_CHECKING("gfortran")
   AC_CHECK_PROG(gfortran, none)

   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-ffree-form -x f95-cpp-input'
   F90FIXED='-ffixed-form -x f95-cpp-input'
   MODFLAG='-M'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static -static-libgfortran'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g ${F90FLAGS}"
       else
	    F90FLAGS="-O3 ${F90FLAGS}"
       fi
   fi

   #do shared specific stuff
   if test "${enable_shared}" = yes ; then
       AC_MSG_CHECKING("rpath based on CXX")
       F90AR=${CXX}
       ARFLAGS='-shared -o'
       AC_DBS_SETUP_RPATH('-Xlinker -rpath', space)
    else
       AR='ar'
    fi

   dnl end of AC_COMPILER_GFORTRAN_F90
])

dnl-------------------------------------------------------------------------dnl
dnl PORTLAND F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_PORTLAND_F90], [dnl

   # Check for working Portland Group F90 compiler

   AC_CHECK_PROG(F90, pgf90, pgf90, none)
   case $F90 in
   *pgf90)
      tmp=`${F90} -V 2>&1 | grep "Portland"`
      if test "${tmp}no" = "no"; then
         AC_MSG_ERROR([not found])
      fi
      ;;
   *) AC_MSG_ERROR([not found]) 
      ;;
   esac

dnl   if test "${F90}" = pgf90 && ${F90} --V 2>&1 | grep "Portland"
dnl   then
dnl       :
dnl   else
dnl       AC_MSG_ERROR([not found])
dnl   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Mfreeform'
   F90FIXED='-Mnofreeform'
   MODFLAG='-module'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -Mbounds -Mchkptr ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_PORTLAND_F90
])

dnl-------------------------------------------------------------------------dnl
dnl COMPAQ F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_COMPAQ_F90], [dnl

   # Check for working compaq F90 compiler

   AC_CHECK_PROG(F90, f95, f95, none)
   if test "${F90}" = f95 && ${F90} -version 2>&1 | grep "Fortran"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE=''
   F90FIXED=''
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-non_shared'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="${F90FREE} -assume byterecl"
       F90FLAGS="${F90FREE} -assume byterecl -automatic -std -warn argument_checking"

       if test "${enable_debug:=no}" = yes
       then
	  # F90FLAGS="-g ${F90FLAGS}"
	    F90FLAGS="-g -check bounds -fpe2 ${F90FLAGS}"
       else
	  # F90FLAGS="-O ${F90FLAGS}"
	    F90FLAGS="-O4 -arch host -assume noaccuracy_sensitive -math_library accurate -tune host ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_COMPAQ_F90
])

dnl-------------------------------------------------------------------------dnl
dnl SUN WORKSHOP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_WORKSHOP_F90], [dnl

   # Check for working WorkShop F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "WorkShop"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # Set F90FREE, F90FIXED, and MODFLAG

   F90FREE='-free'
   F90FIXED='-fixed'
   MODFLAG='-M'

   # Set LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-Bstatic'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_WORKSHOP_F90
])

dnl-------------------------------------------------------------------------dnl
dnl CRAY_F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_CRAY_F90], [dnl

   # Check for working Cray F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # FREE, FIXED AND MODULE FLAGS

   F90FREE='-f free'
   F90FIXED='-f fixed'
   MODFLAG='-p'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	   F90FLAGS="-g ${F90FLAGS}"
       else
	   F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_CRAY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl IRIX MIPS F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_MIPS_F90], [dnl

   # Look for working MIPS compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -version 2>&1 | grep "MIPS"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # Set F90FREE, F90FIXED, and MODFLAG

   F90FREE='-freeform'
   F90FIXED='-col72'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
	#F90FLAGS="${F90FREE} -OPT:Olimit=0"
	F90FLAGS="${F90FREE} -mips4 -r10000 -DEBUG:fullwarn=ON:woff=878,938,1193,1438"

	if test "${enable_debug:=no}" = yes
	then
	  # F90FLAGS="-g ${F90FLAGS}"
	    F90FLAGS="-g -check_bounds -DEBUG:trap_uninitialized=ON ${F90FLAGS}"
	else
	  # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	    F90FLAGS="-O3 -OPT:IEEE_arithmetic=2:roundoff=2 ${F90FLAGS}"
	fi
   fi

   dnl end of AC_COMPILER_MIPS_F90
])

dnl-------------------------------------------------------------------------dnl
dnl HP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_HP_F90], [dnl

   # CHECK FOR WORKING HP F90 COMPILER
   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} +version 2>&1 | grep "HP"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='+source=free'
   F90FIXED='+source=fixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)
   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='+noshared'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} +U77"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -C ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_HP_F90
])

dnl-------------------------------------------------------------------------dnl
dnl INTEL F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_INTEL_F90], [dnl

   # CHECK FOR WORKING INTEL F90 COMPILER
   AC_CHECK_PROG(F90, ifort, ifort, none)
   if test `basename ${F90}` = ifort && ${F90} -V 2>&1 | grep "Intel"; then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='-FR'
   F90FIXED='-FI'
   MODFLAG='-I '
   MODSUFFIX='mod'

   # LINKER AND LIBRARY (AR)
   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} -e95"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -C -implicitnone ${F90FLAGS}"
       else
	    F90FLAGS="-O3 -fno-alias -tpp7 -ipo -pad -align ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_INTEL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl NAG F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_NAG_F90], [dnl

   # CHECK FOR WORKING NAG F90 COMPILER
   AC_CHECK_PROG(F90, f95, f95, none)
   if test "${F90}" = f95 && ${F90} -V 2>&1 | grep "NAGWare"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='-free'
   F90FIXED='-fixed'
   MODFLAG='-I '

   # LINKER AND LIBRARY (AR)
   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-unsharedf95'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} -colour -info -target=native"

       if test "${enable_debug:=no}" = yes
       then
          # only use first line if memory error is suspected, too much output
          #   otherwise
	  # F90FLAGS="-g -C -mtrace=size -nan -u ${F90FLAGS}"
	    F90FLAGS="-g -C -nan -u ${F90FLAGS}"
       else
	    F90FLAGS="-O4 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_NAG_F90
])

dnl-------------------------------------------------------------------------dnl
dnl ABSOFT F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_COMPILER_ABSOFT_F90], [dnl

   # CHECK FOR WORKING ABSOFT F90 COMPILER
   AC_CHECK_PROG(F90, f95, f95, none)
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE=''
   F90FIXED=''
   MODFLAG='-p '

   # LINKER AND LIBRARY (AR)
   F90LD='${F90}'
   F90AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC=''

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="-cpu:host -en"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -et -m0 -M399,1193,878 -Rb -Rc -Rs -Rp -trap=ALL ${F90FLAGS}"
       else
	    F90FLAGS="-O3 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_ABSOFT_F90
])

dnl ========================================================================
