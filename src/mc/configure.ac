dnl-------------------------------------------------------------------------dnl
dnl draco/src/mc configure.in
dnl Thomas M. Evans
dnl Tue Apr 13 17:58:54 1999
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for mc
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([OS_Mesh.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(mc/config.h:config.h.in)
package="mc"

dnl VENDOR SETUPS
AC_MPI_SETUP(pkg)
AC_SPRNG_SETUP(pkg)

dnl installation
AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(viz meshReaders rng c4 traits ds++)

AC_NEEDS_LIBS_TEST(RTT_Format_Reader)

dnl runtesting
AC_RUNTESTS(tstConfig tstMath tstCoord tstLayout tstOSMesh \
            tstFields tstTopology tstParallel_Data_Op \
	    tstTET_Mesh tstSampler tstRZWedge_Mesh \
            tstComm_Patterns tstOSMesh_Pack \
            tstRZWedge_Mesh_Pack tstParticle_Buffer \
            tstCommunicator tstParticle_IO tstConstants, 1 2 4)

dnl DRACO SETUP
AC_DRACO_ENV

dnl do eight byte word setup
AC_DETERMINE_WORD_TYPES
AC_CHECK_EIGHT_BYTE_INT_TYPE

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
                 test/OS_Input:test/OS_Input.in \
	         test/OS_Input_3D:test/OS_Input_3D.in \
                 test/RZWedge_Input:test/RZWedge_Input.in])
AC_OUTPUT

# copy TET meshes to test directory (because they are big)
echo copying ${srcdir}/test/TET_MODEL_1.in to test/TET_MODEL_1
cp ${srcdir}/test/TET_MODEL_1.in test/TET_MODEL_1
echo copying ${srcdir}/test/TET_RTT_1.in to test/TET_RTT_1
cp ${srcdir}/test/TET_RTT_1.in test/TET_RTT_1
echo copying ${srcdir}/test/TET_RTT_2.in to test/TET_RTT_2
cp ${srcdir}/test/TET_RTT_2.in test/TET_RTT_2
echo copying ${srcdir}/test/TET_RTT_3.in to test/TET_RTT_3
cp ${srcdir}/test/TET_RTT_3.in test/TET_RTT_3
echo copying ${srcdir}/test/TET_RTT_4.in to test/TET_RTT_4
cp ${srcdir}/test/TET_RTT_4.in test/TET_RTT_4
echo copying ${srcdir}/test/TET_RTT_5.in to test/TET_RTT_5
cp ${srcdir}/test/TET_RTT_5.in test/TET_RTT_5

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/mc configure.in
dnl-------------------------------------------------------------------------dnl
