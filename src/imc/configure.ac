dnl-------------------------------------------------------------------------dnl
dnl draco/src/imc configure.in
dnl Thomas M. Evans
dnl Mon Jan 11 17:09:11 1999
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl Copyright © 2003 The Regents of the University of California.
dnl-------------------------------------------------------------------------dnl
dnl configure.in for imc
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Source.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(imc/config.h:config.h.in)
package="imc"

# DRACO VENDOR SETUPS
AC_MPI_SETUP(pkg)
AC_SPRNG_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(mc viz meshReaders cdi rng c4 traits ds++)
AC_NEEDS_LIBS_TEST(cdi_analytic)

# testing setup
AC_RUNTESTS(tstFrequency \
            tstMat_State \
            tstOpacity \
            tstOpacity_Builder_Helper \
	    tstFlat_Mat_State_Builder \
            tstCDI_Mat_State_Builder \
            tstMesh_Operations \
            tstSource_Builder \
            tstTally \
            tstRandom_Walk \
            tstAzimuthal_Mesh \
            tstExtrinsic_Tracker \
            tstSurface_Sub_Tally \
            tstSurface_tracker \
            tstRZ_Surface_tracking \
	    tstParticle \
            tstParticle_trans \
            tstTransporter, 1 2 4)

# DRACO SETUP
AC_DRACO_ENV

# OUTPUTS
AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in\
                 test/OS_Input:test/OS_Input.in\
	         test/OS_Input_3D:test/OS_Input_3D.in\
                 test/RZWedge_Input:test/RZWedge_Input.in])
AC_OUTPUT 

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/imc configure.in
dnl-------------------------------------------------------------------------dnl
