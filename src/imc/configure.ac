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
            tstRW_Sub_Tally \
            tstRandom_Walk \
            tstAzimuthal_Mesh \
            tstExtrinsic_Tracker \
            tstSurface_Sub_Tally \
            tstSurface_tracker \
            tstRZ_Surface_tracking \
	    tstParticle \
            tstParticle_trans \
            tstIMC_Objects_Builder \
            tstTransporter, 1 2 4)

# DRACO SETUP
AC_DRACO_ENV

# OUTPUTS
AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])

AC_OUTPUT 

# TEST MESHES
echo copying ${srcdir}/test/OS_Input.in to test/OS_Input
cp ${srcdir}/test/OS_Input.in test/OS_Input
echo copying ${srcdir}/test/OS_Input_3D.in to test/OS_Input_3D
cp ${srcdir}/test/OS_Input_3D.in test/OS_Input_3D
echo copying ${srcdir}/test/OS_Input_B.in to test/OS_Input_B
cp ${srcdir}/test/OS_Input_B.in test/OS_Input_B
echo copying ${srcdir}/test/RZWedge_Input.in to test/RZWedge_Input
cp ${srcdir}/test/RZWedge_Input.in test/RZWedge_Input
echo copying ${srcdir}/test/Sphyramid_Input.in to test/Sphyramid_Input
cp ${srcdir}/test/Sphyramid_Input.in test/Sphyramid_Input

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/imc configure.in
dnl-------------------------------------------------------------------------dnl
