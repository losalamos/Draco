dnl-------------------------------------------------------------------------dnl
dnl draco/src/meshReaders_Services configure.in
dnl Todd Adams
dnl Tue Sep 21 15:09:59 1999
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for meshReaders_Services
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Connect.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(meshReaders_Services/config.h:config.h.in)
package="meshReaders_Services"

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(meshReaders ds++)
AC_NEEDS_LIBS_TEST(RTT_Format_Reader)

AC_RUNTESTS(TestmeshReadersServices, scalar)

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in   \
                 test/Makefile:../../config/Makefile.test.in \
                 test/rttdef.mesh:test/rttdef.mesh.in\
	         test/rttamr.mesh:test/rttamr.mesh.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/meshReaders_Services/configure.in
dnl-------------------------------------------------------------------------dnl
