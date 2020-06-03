//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss test
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "cdi_ndi/NDI_CP_Eloss.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/SystemCall.hh"
#include "ds++/dbc.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using rtt_cdi_ndi::NDI_CP_Eloss;
using rtt_dsxx::soft_equiv;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

void gendir_test(rtt_dsxx::UnitTest &ut) {
  rtt_cdi::CParticle target(1001, 0.);
  rtt_cdi CParticle projectile(2004, 0.);

  // Write a custom gendir file to deal with NDI-required absolute path to data
  NDI_CP_Eloss eloss("gendir path", "dedx", target, projectile);o

  if (ut.numFails == 0) {
    PASSMSG("NDI_CP_Eloss test passes.");
  } else {
    FAILMSG("NDI_CP_Eloss test fails.");
  }
}

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    gendir_test(ut);
    //std::string gendir_default = rtt_dsxx::getFilenameComponent(
    //    std::string(NDI_DATA_DIR) + rtt_dsxx::dirSep + "gendir",
    //    rtt_dsxx::FilenameComponent::FC_NATIVE);

    //if (rtt_dsxx::fileExists(gendir_default)) {
    //  gendir_default_test(ut);
    //}
  }
  UT_EPILOG(ut);
}
