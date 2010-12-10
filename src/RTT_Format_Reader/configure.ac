dnl-------------------------------------------------------------------------dnl
dnl draco/src/RTT_Format_Reader configure.in
dnl Todd Adams
dnl Tue Sep 21 15:09:59 1999
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for RTT_Format_Reader
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([RTT_Format_Reader.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(RTT_Format_Reader/config.h:config.h.in)
package="RTT_Format_Reader"

dnl VENDOR SETUPS
AC_MPI_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(meshReaders \
              mesh_element \
              c4 \
              ds++)
AC_RUNTESTS(TestRTTMeshReader, scalar)
AC_RUNTESTS(TestRTTFormatReader, 1)

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_0_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in   \
                 test/Makefile:../../config/Makefile.test.in \
                 test/rttdef.mesh:test/rttdef.mesh.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/RTT_Format_Reader configure.in
dnl-------------------------------------------------------------------------dnl
