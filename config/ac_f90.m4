dnl ========================================================================
dnl 
dnl 	Author:	Mark G. Gray
dnl 		Los Alamos National Laboratory
dnl 	Date:	Sun Apr  2 12:33:41 MDT 2000
dnl 
dnl 	Copyright (c) 2000 U. S. Department of Energy. All rights reserved.
dnl 
dnl	$Id$
dnl 
dnl ========================================================================

dnl NAME

dnl	AC_LANG_F90, AC_PROG_GM4, AC_REQUIRE_GM4, AC_PROG_F90 - 
dnl     Fortran 90 and gm4 macros for autoconf 

dnl SYNOPSIS/USAGE

dnl	AC_LANG_F90
dnl 	AC_PROG_GM4
dnl 	AC_REQUIRE_GM4
dnl	AC_PROG_F90

dnl DESCRIPTION

dnl	AC_LANG_F90 sets up the compile and link test invocation using
dnl     F90 and the extension F90EXT.

dnl 	AC_PROG_GM4 sets the output variable GM4 to a command that runs 
dnl	the GNU m4 preprocessor.  Looks for gm4 and then m4, and verifies 
dnl	the gnu version by running ${GM4} --version

dnl	AC_REQUIRE_GM4 ensures that GM4 has been found

dnl	AC_PROG_F90 determines the Fortran 90 compiler to use and the
dnl     appropriate extension for free format source code.  If F90
dnl	is not already set in the environment, check for `f90'.  Try
dnl     the compiler with the extensions `f90', `F', `F90', and `f'.
dnl     Set the output variable F90 to the name of the compiler found
dnl     and F90EXT to the name of the free format source extension.
dnl     If the output variable F90FLAGS was not already set, set it to
dnl      `-g'.  If the output variables F90FREE and F90FIXED were not
dnl	already set, try to guess their values from the target.  If
dnl	the output variable MODFLAG was not already set, try to guess
dnl	its value from the target.

dnl BUGS

dnl	These macros have only been tested on a limited number of
dnl	machines.   AC_PROG_F90 can fail due to vendor non-standard
dnl	file extentions or incorrect free/fixed source defaults.
dnl	F90FREE and F90FIXED correctly set for only a few known
dnl	targets.  AC_F90_MOD can be confused by other files created
dnl	during compilation.  As with other autoconf macros, any file
dnl	named [Cc]onftest* will be overwritten!

dnl See acgeneral.m4...

dnl ### Define output for F90 compiler mesages

define(AC_FD_F90, 5)dnl  See AC_FD_CC
[#] AC_FD_F90 compiler messages saved in config.log

dnl Restore the current language from the stack.
dnl AC_LANG_RESTORE()
popdef([AC_LANG_RESTORE])dnl remove acgeneral.m4 definition
pushdef([AC_LANG_RESTORE], [dnl See AC_LANG_RESTORE
ifelse(AC_LANG_STACK, [C], [AC_LANG_C],dnl
AC_LANG_STACK, [CPLUSPLUS], [AC_LANG_CPLUSPLUS],dnl
AC_LANG_STACK, [F90], [AC_LANG_F90],dnl add F90 to list
AC_LANG_STACK, [FORTRAN77], [AC_LANG_FORTRAN77])[]popdef([AC_LANG_STACK])])

dnl ### Selecting which language to use for testing
AC_DEFUN(AC_LANG_F90, [dnl See AC_LANG_CC
   define([AC_LANG], [F90])dnl
   ac_compile='${F90-f90} -c $F90FLAGS $F90FREE Conftest.$F90EXT 1>&AC_FD_F90'
   ac_link='${F90-f90} -o Conftest $F90FLAGS $F90FREE $LDFLAGS Conftest.$F90EXT $LIBS 1>&AC_FD_F90'
   cross_compiling=$ac_cv_prog_f90_cross
])

dnl See acspecific.m4...

AC_DEFUN(AC_PROG_GM4, [dnl
   AC_MSG_CHECKING(how to run the Gnu m4 preprocessor)
   AC_CHECK_PROGS(GM4, gm4 m4, none)
   if test "${GM4}" != none && ${GM4} --version 2>&1 | grep "GNU"
   then
       AC_MSG_RESULT([found])
   else
       AC_MSG_ERROR([not found])
   fi
])

dnl Require finding the gnu m4 preprocessor if F90 is the current language
AC_DEFUN(AC_REQUIRE_GM4, [dnl See AC_REQUIRE_CPP
   ifelse(AC_LANG, F90, [AC_REQUIRE([AC_PROG_GM4])])
]) 

dnl ### Checks for programs

AC_DEFUN(AC_PROG_F90, [dnl
   AC_REQUIRE([AC_LANG_F90])
   AC_REQUIRE([AC_DRACO_ENV])
   AC_PROG_F90_WORKS
   AC_F90_MOD
])

AC_DEFUN(AC_PROG_F90_WORKS, [dnl
   AC_MSG_CHECKING([whether the F90 compiler ($F90 $F90FLAGS $F90FREE $LDFLAGS) works])
   AC_LANG_SAVE
   AC_LANG_F90
   for F90EXT in f90 F F90 f
   do
      AC_TRY_F90_COMPILER([program main; end], 
                    ac_cv_prog_f90_works, ac_cv_prog_f90_cross)
      if test $ac_cv_prog_f90_works = yes; then
          break
      fi
   done
   AC_LANG_RESTORE
   AC_MSG_RESULT($ac_cv_prog_f90_works)
   if test $ac_cv_prog_f90_works = no; then
      AC_MSG_ERROR([installation or configuration problem: F90 compiler cannot create executables.])
   fi
   AC_MSG_CHECKING([whether the F90 compiler ($F90 $F90FLAGS $F90FREE $LDFLAGS) is a cross-compiler])
   AC_MSG_RESULT($ac_cv_prog_f90_cross)
   cross_compiling=$ac_cv_prog_f90_cross
])

AC_DEFUN(AC_TRY_F90_COMPILER, [dnl
   cat > Conftest.$F90EXT <<EOF
[$1]
EOF
   if AC_TRY_EVAL(ac_link) && test -s Conftest; then
      [$2]=yes
      # If we can't run a trivial program, 
      # we are probably using a cross compiler.
      if (./Conftest; exit) 2>/dev/null; then
         [$3]=no
      else
         [$3]=yes
      fi
   else
      echo "configure: failed program was:" >&AC_FD_CC
      cat Conftest.$F90EXT >&AC_FD_CC
      [$2]=no
   fi
   rm -fr Conftest*
   AC_SUBST(F90EXT)dnl
   AC_SUBST(F90FLAGS)dnl
   AC_SUBST(F90FREE)dnl
   AC_SUBST(F90FIXED)dnl
])

dnl ### Checks for module information

AC_DEFUN(AC_F90_MOD,[dnl
   AC_MSG_CHECKING([the F90 compiler module name])
   AC_LANG_SAVE
   AC_LANG_F90
   if test "" = "${MODNAME}" -a "" = "${MODSUFFIX}"; then
      rm -f conftest*
      cat > Conftest.${F90EXT} <<EOF
module Conftest_foo
end module Conftest_foo
EOF
   if AC_TRY_EVAL(ac_compile) && test -s Conftest.o; then
      rm -f Conftest.${F90EXT} Conftest.o
      modfile=`ls | grep -i Conftest`
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
      *)             MODSUFFIX=o; MODNAME=Filename ;;
      esac
      else
         echo "configure: failed program was:" >&AC_FD_F90
         cat Conftest.${F90EXT} >&AC_FD_F90
      fi
      rm -f $modfile
   fi
   AC_MSG_RESULT($MODNAME.$MODSUFFIX)
   AC_SUBST(MODSUFFIX)dnl
   AC_SUBST(MODNAME)dnl
   AC_SUBST(MODFLAG)dnl
   AC_LANG_RESTORE
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_f90.m4
dnl-------------------------------------------------------------------------dnl
