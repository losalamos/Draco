//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   nuclear_data/test/tNuclear_Data.cc
 * \author B. R. Ryan
 * \date   2019 Jun 24
 * \brief  Implementation file for tNuclear_Data
 * \note   Copyright (C) 2017-2019 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "nuclear_data/Nuclear_Data.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
//#include <sstream>
#include <stdexcept>

namespace rtt_nuclear_data_test {

using rtt_dsxx::soft_equiv;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

void nuclear_data_file_test(rtt_dsxx::UnitTest &ut) {

  std::cout << "\n---------------------------------------------------------\n"
            << "   Test Draco code calling Nuclear_Data routines\n"
            << "---------------------------------------------------------\n";

  // Open a sample photoatomic ACE datafile lifted from ACEtk
  const std::string filename = ut.getTestSourcePath() + "1000.14p";
  std::cout << "Attempting to construct a Nuclear_Data object...\n"
            << std::endl;
  std::unique_ptr<rtt_nuclear_data::Nuclear_Data> nuclear_data;

  try {
    nuclear_data.reset(new rtt_nuclear_data::Nuclear_Data(filename));
  } catch (...) {
    FAILMSG("Failed to construct a Nuclear_Data object!");
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Check some of the parsed data
  const std::vector<double> ESZG = nuclear_data->get_ESZG();
  const std::vector<double> LNEPS = nuclear_data->get_LNEPS();
  const std::vector<double> SWD = nuclear_data->get_SWD();

  if (!soft_equiv(ESZG[0], -13.75744544419))
    FAILMSG("ESZG[0] read incorrectly!");
  if (!soft_equiv(ESZG[3234], 0.))
    FAILMSG("ESZG[3234] read incorrectly!");

  if (!soft_equiv(LNEPS[0], 1.4e-05))
    FAILMSG("LNEPS[0] read incorrectly!");

  if (!soft_equiv(SWD[0], 31.))
    FAILMSG("SWD[0] read incorrectly!");
  if (!soft_equiv(SWD[8401], 12987.1))
    FAILMSG("SWD[8401] read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read ESZG, LNEPS, and SWD!" << std::endl;
  }

  if (ut.numFails == 0) {
    PASSMSG("Successfully linked Draco against Nuclear_Data.");
  } else {
    FAILMSG("Did not successfully link Draco against Nuclear_Data.");
  }
}

void nuclear_data_bad_file_test(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n---------------------------------------------------------\n"
            << "    Test Nuclear_Data bad file handling    \n"
            << "---------------------------------------------------------\n";

  std::string filename = ut.getTestSourcePath() + "nonexistent_ace_file";
  std::cout << "Testing with a non-existent file...\n" << std::endl;
  std::unique_ptr<rtt_nuclear_data::Nuclear_Data> nuclear_data;

  bool caught = false;
  try {
    nuclear_data.reset(new rtt_nuclear_data::Nuclear_Data(filename));
  } catch (rtt_dsxx::assertion &asrt) {
    std::cout << "Draco exception thrown: " << asrt.what() << std::endl;
    caught = true;
  } catch (const int &) {
    std::cout << "Nuclear data exception thrown. " << std::endl;
    caught = true;
  }

  if (!caught)
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully caught a Nuclear_Data constructor exception.");
  } else {
    FAILMSG("Did not successfully catch a Nuclear_Data constructor exception.");
  }
}

} // namespace rtt_nuclear_data_test

//----------------------------------------------------------------------------//
// Main
//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    rtt_nuclear_data_test::nuclear_data_file_test(ut);
    rtt_nuclear_data_test::nuclear_data_bad_file_test(ut);
  }
  UT_EPILOG(ut);

  return 0;

  try {
    rtt_nuclear_data_test::nuclear_data_file_test(ut);
  } catch (rtt_dsxx::assertion &err) {
    std::string msg = err.what();
    if (msg != std::string("Success")) {
      std::cout << "ERROR: While testing " << argv[0] << ", " << err.what()
                << std::endl;
      return 1;
    }
    return 0;
    /*} catch (rtt_dsxx::exception &err) {
    std::cout << "ERROR: While testing " << argv[0] << ", " << err.what() 
              << std::endl;
    return 1;*/
  } catch (...) {
    std::cout << "ERROR: While testing " << argv[0] << ", "
              << "An unknown exception was thrown" << std::endl;
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------//
// end of tNuclear_Data.cc
//----------------------------------------------------------------------------//
