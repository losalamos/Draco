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
dnl     `f95', and `xlf90', in that order.  Set the output variable `F90' 
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

dnl See acgeneral.m4...

dnl Restore the current language from the stack.

dnl AC_LANG_RESTORE()
popdef([AC_LANG_RESTORE])dnl remove acgeneral.m4 definition
pushdef([AC_LANG_RESTORE], [dnl See AC_LANG_RESTORE
   indir([AC_LANG_]AC_LANG_STACK)
   popdef([AC_LANG_STACK])
])

dnl ### Selecting which language to use for testing
dnl     See AC_LANG_C, AC_LANG_CPLUSPLUS, AC_LANG_FORTRAN77
AC_DEFUN(AC_LANG_FORTRAN90, [dnl 
   define([AC_LANG], [FORTRAN90])dnl
   ac_ext=f90
   ac_compile='${F90-f90} -c $F90FLAGS conftest.$ac_ext 1>&AC_FD_CC'
   ac_link='${F90-f90} -o conftest${ac_exeext} $F90FLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&AC_FD_CC'
   cross_compiling=$ac_cv_prog_f90_cross
])

dnl See acspecific.m4...

dnl ### Checks for programs

AC_DEFUN(AC_PROG_F90, [dnl
   if test -z "$F90"
   then
       AC_CHECK_PROGS(F90, f90 F90 f95 xlf90)
       test -z "$F90" && AC_MSG_ERROR([no acceptable Fortran 90 compiler found in \$PATH])
   fi
   AC_PROG_F90_WORKS
   AC_PROG_F90_MOD

   dnl Check whether -g works, even if F90FLAGS is set, in case the package
   dnl plays around with F90FLAGS (such as to build both debugging and
   dnl normal versions of a library), tasteless as that idea is.
   ac_test_F90FLAGS="${F90FLAGS+set}"
   ac_save_F90FLAGS="$F90FLAGS"
   AC_PROG_F90_G
   if test "$ac_test_F90FLAGS" = set 
   then
       F90FLAGS="$ac_save_F90FLAGS"
   elif test $ac_cv_prog_f90_g = yes 
   then
       F90FLAGS="-g"
   else
       F90FLAGS=
   fi
])

AC_DEFUN(AC_PROG_F90_WORKS, [dnl
   AC_MSG_CHECKING([whether the F90 compiler ($F90 $F90FLAGS $LDFLAGS) works])
   AC_LANG_SAVE
   AC_LANG_FORTRAN90
   AC_TRY_F90_COMPILER([program main; end], 
	ac_cv_prog_f90_works, ac_cv_prog_f90_cross)
   AC_LANG_RESTORE
   AC_MSG_RESULT($ac_cv_prog_f90_works)
   if test $ac_cv_prog_f90_works = no; then
      AC_MSG_ERROR([installation or configuration problem: F90 compiler cannot create executables.])
   fi
   AC_MSG_CHECKING([whether the F90 compiler ($F90 $F90FLAGS $LDFLAGS) is a cross-compiler])
   AC_MSG_RESULT($ac_cv_prog_f90_cross)
   cross_compiling=$ac_cv_prog_f90_cross
])

dnl Test whether the Fortran 90 compiler can accept the `-g' option to
dnl enable debugging.
dnl 
dnl AC_PROG_F90_G()
AC_DEFUN(AC_PROG_F90_G, [dnl
   AC_LANG_SAVE
   AC_LANG_FORTRAN90
   AC_CACHE_CHECK(whether $F90 accepts -g, ac_cv_prog_f90_g, [dnl
      cat > conftest.$ac_ext << EOF
subroutine conftest()
end
EOF
      if test -z "`$F90 -g -c conftest.$ac_ext 2>&1`"; then
          ac_cv_prog_f90_g=yes
      else
          ac_cv_prog_f90_g=no
      fi
      rm -f conftest*
   ])
   AC_LANG_RESTORE
])

dnl Try to compile, link and execute TEST-PROGRAM.  Set WORKING-VAR to
dnl `yes' if the current compiler works, otherwise set it ti `no'.  Set
dnl CROSS-VAR to `yes' if the compiler and linker produce non-native
dnl executables, otherwise set it to `no'.  Before calling
dnl `AC_TRY_COMPILER()', call `AC_LANG_*' to set-up for the right
dnl language.
dnl 
dnl AC_TRY_F90_COMPILER(TEST-PROGRAM, WORKING-VAR, CROSS-VAR)
AC_DEFUN(AC_TRY_F90_COMPILER,
[cat > conftest.$ac_ext << EOF
[$1]
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  [$2]=yes
  # If we can't run a trivial program, we are probably using a cross compiler.
  if (./conftest${ac_exeext}; exit) 2>/dev/null; then
    [$3]=no
  else
    [$3]=yes
  fi
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
  [$2]=no
fi
rm -fr conftest*])


dnl ### Checks for module information

AC_DEFUN(AC_PROG_F90_MOD,[dnl
   AC_MSG_CHECKING([the F90 compiler module name])
   if test -z "$MODNAME" -a -z "$MODSUFFIX"
   then
       AC_LANG_SAVE
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
       AC_LANG_RESTORE
   fi
   AC_MSG_RESULT($MODNAME.$MODSUFFIX)
   AC_SUBST(MODSUFFIX)dnl
   AC_SUBST(MODNAME)dnl
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_f90.m4
dnl-------------------------------------------------------------------------dnl
