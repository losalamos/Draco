dnl-------------------------------------------------------------------------dnl
dnl draco/src/ds++ configure.in
dnl Thomas M. Evans
dnl Thu Jan  7 16:58:31 1999
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for ds++
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([SP.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(ds++/config.h:config.h.in)
package="ds++"

# Toggle between autoconf and cmake build systems when constructing
# config.h from config.h.in
AC_DEFINE(AUTOCONF autoconf)

AC_INSTALL_HEADERS
AC_INSTALL_LIB

AC_RUNTESTS(tstAllocators \
            tstArray \
            tstAssert \     
            tstBounds \
            tstCheck_Strings \
            tstData_Table \
            tstSafe_Ptr \
            tstDBC_Array \
            tstDBC_Ptr \
            tstDynArray \
            tstEndian \
            tstFile_Streams \
            tstIndex_Converter \
            tstIndex_Counter \
            tstIndex_Set \
            tstMat \
            tstMat1RA \
            tstMat2RA \
            tstMat3RA \
            tstMat4RA \
            tstMat5RA \
            tstPacking_Utils \
            tstRCF \
            tstRange_finder \
            tstSP \
            tstSafe_Divide \
            tstScalarUnitTest \
            tstSoft_Equiv \
            tstSortPermutation \
            tstVector_Lite \
            tstSlice \
            tstisFinite \
            tstsquare \
            tstpythag \
            tstabs \
            tstconj \
            tstcube \
            tstsign \
            tstvalue \
            tstField_Traits \
            tstto_string \
            tstdbc \
            tstHomogeneous_New \
            , scalar)

AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-5_21_0)

dnl -------------------------------------------
dnl look for flexelint, dot and ps2pdf on Linux
dnl -------------------------------------------
AC_PATH_PROG(FLEXELINT_PATH, lin, null, /usr/local/flexelint/bin )
AC_PATH_PROG(DOTCMD, dot, null)
AC_PATH_PROG(PS2PDF, ps2pdf, null)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
                 doc/Makefile:../../config/Makefile.doc.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/ds++ configure.in
dnl-------------------------------------------------------------------------dnl
