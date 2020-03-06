//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/test/tstNDI_AMW.cc
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AMW test
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "cdi/CDI.hh"
#include "cdi_ndi/NDI_AMW.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using rtt_cdi_ndi::NDI_AMW;
using rtt_dsxx::soft_equiv;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

void amw_test(rtt_dsxx::UnitTest &ut) {

  // Write a custom gendir file to deal with NDI-required absolute path to data
  std::string gendir_in = "gendir_tmp.all";
  std::string gendir_tmp_path = ut.getTestInputPath() + gendir_in;
  std::string data_path = ut.getTestSourcePath() + "ndi_data";
  std::ofstream gendir_tmp_file;
  gendir_tmp_file.open(gendir_tmp_path);
  gendir_tmp_file << "multigroup_neutron\n";
  gendir_tmp_file << "  z=1001.700nm  d=11/19/2010  l=mendf71x\n";
  gendir_tmp_file << "    f=fake/data/path\n";
  gendir_tmp_file << "    ft=asc  ln=2  o=28\n";
  gendir_tmp_file << "    ng=618  t=2.5300642359999999e-08  s0=10000000000\n";
  gendir_tmp_file
      << "    aw=1.0078249887344399  awr=0.99916729999999998  end\n";
  gendir_tmp_file.close();

  NDI_AMW ndi_amw(gendir_tmp_path);

  int proton_zaid = 1001;
  double proton_amw = ndi_amw.get_amw(proton_zaid);

  FAIL_IF_NOT(soft_equiv(proton_amw, 1.0078249887344399, 1.e-8));

  if (ut.numFails == 0) {
    PASSMSG("NDI_AMW test passes.");
  } else {
    FAILMSG("NDI_AMW test fails.");
  }
}

//----------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    amw_test(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// end of tstNDI_AMW.cc
//----------------------------------------------------------------------------//
