dnl ========================================================================
dnl 
dnl 	Author:	Mark G. Gray
dnl 		Los Alamos National Laboratory
dnl 	Date:	Sun Apr  2 12:33:41 MDT 2000
dnl 
dnl 	Copyright (c) 2000 Free Software Foundation
dnl 
dnl	$Id$
dnl 
dnl ========================================================================

dnl NAME

dnl	AC_LANG_FORTRAN90, AC_PROG_F90, 
dnl     Fortran 90 macros for autoconf 

dnl SYNOPSIS/USAGE

dnl	AC_LANG_FORTRAN90
dnl	AC_PROG_F90

dnl DESCRIPTION

dnl	AC_LANG_FORTRAN90 sets up the compile and link test.  Use F90, 
dnl     F90FLAGS, and LDFLAGS for test programs.

dnl	AC_PROG_F90 determines a Fortran 90 compiler to use.  If F90
dnl	is not already set in the environment, check for `f90', `F90',
dnl     `f95', and `xlf95', in that order.  Set the output variable `F90' 
dnl     to the name of the compiler found. 

dnl     If the output variable F90FLAGS was not already set, set it to
dnl      `-g'.  
dnl
dnl     If MODNAME and MODSUFFIX are not already set in the environment, 
dnl     test for MODNAME and MODSUFFIX.  Set the output variables `MODNAME'
dnl     and `MODSUFFIX' to the module name and suffix conventions, 
dnl     respectively.

dnl BUGS

dnl	These macros have only been tested on a limited number of
dnl	machines.   AC_PROG_F90 can fail due to vendor non-standard
dnl	file extentions or incorrect free/fixed source defaults.
dnl	F90FREE and F90FIXED correctly set for only a few known
dnl	targets.  AC_F90_MOD can be confused by other files created
dnl	during compilation.  As with other autoconf macros, any file
dnl	named [Cc]onftest* will be overwritten!

dnl ### Selecting which language to use for testing
dnl     See AC_LANG_C, AC_LANG_CPLUSPLUS, AC_LANG_FORTRAN77
AC_DEFUN(AC_LANG_FORTRAN90, [dnl 
   define([AC_LANG], [FORTRAN90])dnl
   ac_ext=f90
   ac_compile='${F90-f90} -c $F90FLAGS conftest.$ac_ext 1>&AC_FD_CC'
   ac_link='${F90-f90} -o conftest${ac_exeext} $F90FLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&AC_FD_CC'
   cross_compiling=$ac_cv_prog_f90_cross
])


dnl ### Checks for module information

AC_DEFUN(AC_PROG_F90_MOD,[dnl
   AC_MSG_CHECKING([the F90 compiler module name])
   if test -z "$MODNAME" -a -z "$MODSUFFIX"
   then
       AC_LANG_FORTRAN90
       rm -f conftest*
       cat > conftest.$ac_ext <<EOF
module conftest_foo
end module conftest_foo
EOF
       if AC_TRY_EVAL(ac_compile) && test -s conftest.o 
       then
           rm -f conftest.$ac_ext conftest.o
           modfile=`ls | grep -i conftest`
           test "${modfile+set}" = set || AC_MSG_ERROR([unknown modfile: set MODSUFFIX and MODNAME in environment])
           MODSUFFIX=`expr "$modfile" : ".*\.\(.*\)"`
           MODNAME=`basename $modfile .$MODSUFFIX`
           case "$MODNAME" in
           conftest)      MODNAME=filename ;;
           Conftest)      MODNAME=Filename ;;
           CONFTEST)      MODNAME=FILENAME ;;
           conftest_foo)  MODNAME=modname ;;
           Conftest_foo)  MODNAME=Modname ;;
           CONFTEST_FOO)  MODNAME=MODNAME ;;
           *)             MODNAME=Filename 
                          MODSUFFIX=o ;;
           esac
       else
           echo "configure: failed program was:" >&AC_FD_CC
           cat conftest.$ac_ext >&AC_FD_CC
       fi
       rm -f $modfile
   fi
   AC_MSG_RESULT($MODNAME.$MODSUFFIX)
   AC_SUBST(MODSUFFIX)dnl
   AC_SUBST(MODNAME)dnl
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_f90.m4
dnl-------------------------------------------------------------------------dnl
