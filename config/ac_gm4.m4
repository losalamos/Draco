dnl  ========================================================================
dnl 
dnl 	Author:	Mark G. Gray
dnl 		Los Alamos National Laboratory
dnl 	Date:	Fri Apr 21 14:06:27 MDT 2000
dnl 
dnl 	Copyright (c) 2000 U. S. Department of Energy. All rights reserved.
dnl 
dnl  ========================================================================

dnl NAME

dnl     AC_PROG_GM4, AC_REQUIRE_GM4 - 
dnl     gm4 macros for autoconf 

dnl SYNOPSIS/USAGE

dnl 	AC_PROG_GM4
dnl 	AC_REQUIRE_GM4

dnl DESCRIPTION

dnl 	AC_PROG_GM4 sets the output variable GM4 to a command that runs 
dnl	the GNU m4 preprocessor.  Looks for gm4 and then m4, and verifies 
dnl	the gnu version by running ${GM4} --version

dnl	AC_REQUIRE_GM4 ensures that GM4 has been found

dnl  ========================================================================
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

AC_DEFUN([AC_PROG_GM4], [dnl
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
AC_DEFUN([AC_REQUIRE_GM4], [dnl See AC_REQUIRE_CPP
   ifelse(AC_LANG, F90, [AC_REQUIRE([AC_PROG_GM4])])
]) 

dnl  ========================================================================
