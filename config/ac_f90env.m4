dnl ========================================================================
dnl 
dnl 	Author:	Mark G. Gray
dnl 		Los Alamos National Laboratory
dnl 	Date:	Wed Apr 19 16:39:19 MDT 2000
dnl 
dnl 	Copyright (c) 2000 U. S. Department of Energy. All rights reserved.
dnl 
dnl	$Id$
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

dnl ### Ensure with_f90 set
AC_DEFUN(AC_WITH_F90, [dnl
   : ${with_f90:=yes}
])

dnl
dnl CHOOSE A F90 COMPILER
dnl

dnl defines --with-f90
AC_ARG_WITH(f90,[dnl
  --with-f90[=XL,WorkShop,Fujitsu,Absoft,Cray,MIPS,Compaq]     choose an F90 compiler
])

AC_DEFUN(AC_F90_ENV, [dnl
   AC_REQUIRE([AC_CANONICAL_SYSTEM])

   case "${with_f90:=yes}" in
   XL)
       AC_COMPILER_XL_F90
   ;;
   Fujitsu)
       AC_COMPILER_FUJITSU_F90
   ;;
   Lahey)
       AC_COMPILER_LAHEY_F90
   ;;
   Portland)
       AC_COMPILER_PORTLAND_F90
   ;;
   Absoft)
	AC_COMPILER_ABSOFT_F90
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
   yes)				# guess compiler from target platform
       case "${target}" in   
       rs6000-ibm-aix*)
           AC_COMPILER_XL_F90
       ;;
       sparc-sun-solaris2.*)
           AC_COMPILER_WORKSHOP_F90
       ;;
       i?86-pc-linux*)
           AC_COMPILER_FUJITSU_F90
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
       alphaev67-dec*)
          AC_COMPILER_COMPAQ_F90
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

   AC_SUBST(F90FREE)
   AC_SUBST(F90FIXED)
   AC_SUBST(F90FLAGS)
   AC_SUBST(MODFLAG)
])

dnl-------------------------------------------------------------------------dnl
dnl IBM XLF90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_XL_F90, [dnl

   # Check for working XL F90 compiler

   AC_CHECK_PROG(F90, xlf90, xlf90, none)
   if test "${F90}" != xlf90
   then
       AC_MSG_ERROR([not found])
   fi
  
   # FREE, FIXED AND MODULE FLAGS

   F90FREE='-qfree=f90'
   F90FIXED='-qfixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
       F90FLAGS="-qextchk -qhalt=s -qarch=pwr2 -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	   trapflags="-qinitauto=FF"
	   trapflags="${trapflags} -qflttrap=overflow:underflow:zerodivide:invalid:enable"
	   trapflags="${trapflags} -qsigtrap"
	   F90FLAGS="-g -d -C ${trapflags} -bloadmap:loadmap.dat ${F90FLAGS}"
       else
	   F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_XL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl FUJITSU F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_FUJITSU_F90, [dnl

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

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # COMPILATION FLAGS

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

AC_DEFUN(AC_COMPILER_LAHEY_F90, [dnl

   # Check for working Lahey F90 compiler

   AC_CHECK_PROG(F90, lf95, lf95, none)
   if test "${F90}" = lf95 && ${F90} --version 2>&1 | grep "Lahey"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Free'
   F90FIXED='-Fixed'
   MODFLAG='-M'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
       F90FLAGS="--f95 ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g --chk(aesux) --chkglobal ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_LAHEY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl PORTLAND F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_PORTLAND_F90, [dnl

   # Check for working Portland Group F90 compiler

   AC_CHECK_PROG(F90, pgf90, pgf90, none)
   if test "${F90}" = pgf90 && ${F90} --V 2>&1 | grep "Portland"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Mfreeform'
   F90FIXED='-Mnofreeform'
   MODFLAG='-module'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC=

   # COMPILATION FLAGS

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
dnl ABSOFT F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_ABSOFT_F90, [dnl

   # Check for working Absoft F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)

   # Stupid Absoft Fortran does not have a version header
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE=''
   F90FIXED=''
   MODFLAG='-p'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
       F90FLAGS=""

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_ABSOFT_F90
])

dnl-------------------------------------------------------------------------dnl
dnl COMPAQ F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_COMPAQ_F90, [dnl

   # Check for working compaq F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -version 2>&1 | grep "Compaq"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-free'
   F90FIXED='-fixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-non_shared'

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g ${F90FLAGS}"
       else
	    F90FLAGS="-O ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_COMPAQ_F90
])

dnl-------------------------------------------------------------------------dnl
dnl SUN WORKSHOP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_WORKSHOP_F90, [dnl

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

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # Set COMPILATION FLAGS

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

AC_DEFUN(AC_COMPILER_CRAY_F90, [dnl

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

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # COMPILATION FLAGS

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

AC_DEFUN(AC_COMPILER_MIPS_F90, [dnl

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

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
	F90FLAGS="${F90FREE} -OPT:Olimit=0"

	if test "${enable_debug:=no}" = yes
	then
	    F90FLAGS="-g ${F90FLAGS}"
	else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	fi
   fi

   dnl end of AC_COMPILER_MIPS_F90
])

dnl ========================================================================




