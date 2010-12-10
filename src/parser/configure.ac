dnl-------------------------------------------------------------------------dnl
dnl draco/src/parser configure.ac
dnl Kent G. Budge
dnl Fri Jan 24 11:13:40 2003
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.ac for parser
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Parse_Table.cc])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(parser/config.h:config.h.in)
package="parser"

dnl VENDOR SETUPS
AC_MPI_SETUP(pkg)

dnl Draco version
AC_DRACO_PREREQ

dnl
dnl What to install
dnl
AC_INSTALL_HEADERS
AC_INSTALL_LIB

dnl 
dnl Dependencies
dnl
AC_NEEDS_LIBS( c4 units \
               ds++ )

dnl
dnl Unit tests
dnl
AC_RUNTESTS(tstToken \
            tstToken_Equivalence \
            tstutilities \
            tstUnit \
            tstString_Token_Stream \
            tstParse_Table \
            tstClass_Parser \
            tstAbstract_Class_Parser \
            tstExpression \
            , scalar)

AC_RUNTESTS(tstFile_Token_Stream \
	    tstParallel_File_Token_Stream \
            driver4tstConsole_Token_Stream \
            , 1 2)

dnl
dnl DRACO SETUP
dnl
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_0_0)

dnl 
dnl Create build system files in pwd.
dnl
AC_CONFIG_FILES(Makefile:../../config/Makefile.package.in \
                test/Makefile:../../config/Makefile.test.in)

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/parser configure.ac
dnl-------------------------------------------------------------------------dnl
