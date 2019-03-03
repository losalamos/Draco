//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/test/ReadOdfIpcressFile.cc
 * \author Seth R. Johnson
 * \date   Thu July 10 2008
 * \brief  Regression test based on odfregression10.ipcress, also checks
 *         packing and unpacking.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "cdi_ipcress_test.hh"
#include "cdi/OpacityCommon.hh"
#include "cdi_ipcress/IpcressFile.hh"
#include "cdi_ipcress/IpcressMultigroupOpacity.hh"
#include "cdi_ipcress/IpcressOdfmgOpacity.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cstdio>

using rtt_cdi_ipcress::IpcressFile;
using rtt_cdi_ipcress::IpcressMultigroupOpacity;
using rtt_cdi_ipcress::IpcressOdfmgOpacity;
using rtt_dsxx::soft_equiv;

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::istringstream;
using std::ostringstream;
using std::string;

typedef std::shared_ptr<IpcressOdfmgOpacity const> SP_Goo;
typedef std::vector<double> vec_d;
typedef std::vector<std::vector<double>> vec2_d;

//---------------------------------------------------------------------------//

namespace benchmarkData {
rtt_cdi::Model const model = rtt_cdi::ROSSELAND;
rtt_cdi::Reaction const reaction = rtt_cdi::ABSORPTION;
int const matID = 19000;

int const numGroups = 10;
double const groupBoundaries[numGroups + 1] = {
    1.00000000e-02, 2.51188643e-02, 6.30957344e-02, 1.58489319e-01,
    3.98107171e-01, 1.00000000e+00, 2.51188643e+00, 6.30957344e+00,
    1.58489319e+01, 3.98107171e+01, 1.00000000e+02};

int const numBands = 8;
double const bandBoundaries[numBands + 1] = {0,    0.03, 0.13, 0.33, 0.63,
                                             0.83, 0.93, 0.98, 1};

double const temp = 0.1;
double const dens = 0.1;

double const opacities[numGroups][numBands] = {
    {
        16775.49059752706, // group 1 band 1
        18182.67295203084, // group 1 band 2
        22057.70880057825, // group 1 band 3
        32355.43222366086, // group 1 band 4
        53310.33190227170, // group 1 band 5
        79445.68918112732, // group 1 band 6
        101395.9876836948, // group 1 band 7
        114753.1693005596  // group 1 band 8
    },
    {
        4850.909447480577, // group 2 band 1
        5129.120853712516, // group 2 band 2
        6120.937372844221, // group 2 band 3
        7545.099763678086, // group 2 band 4
        10129.92201796587, // group 2 band 5
        12247.96826888786, // group 2 band 6
        14456.70409922315, // group 2 band 7
        16124.80485208393  // group 2 band 8
    },
    {
        1373.913002095310, // group 3 band 1
        1458.970012669071, // group 3 band 2
        1527.311165164504, // group 3 band 3
        1828.236115918924, // group 3 band 4
        2970.931473569568, // group 3 band 5
        4995.184609533723, // group 3 band 6
        5589.267483055316, // group 3 band 7
        6427.356960341955  // group 3 band 8
    },
    {
        1504.744680800113, // group 4 band 1
        1558.513082092346, // group 4 band 2
        1655.854703874608, // group 4 band 3
        1955.053368870200, // group 4 band 4
        3035.496933071605, // group 4 band 5
        5622.559019580522, // group 4 band 6
        8602.980940325369, // group 4 band 7
        14806.98815916487  // group 4 band 8
    },
    {
        524.9004052253512, // group 5 band 1
        563.1730976990697, // group 5 band 2
        726.6489469107313, // group 5 band 3
        1182.320779000210, // group 5 band 4
        2939.436823682617, // group 5 band 5
        8280.694060920965, // group 5 band 6
        26036.43434148943, // group 5 band 7
        64160.23655271116  // group 5 band 8
    },
    {
        864.1908520112063, // group 6 band 1
        975.8988646476864, // group 6 band 2
        1335.879237840064, // group 6 band 3
        2300.468927485172, // group 6 band 4
        3331.749691177460, // group 6 band 5
        4187.202472888080, // group 6 band 6
        5098.839348085739, // group 6 band 7
        6576.082956860007  // group 6 band 8
    },
    {
        69.82385027300423, // group 7 band 1
        81.26031848373785, // group 7 band 2
        118.3339357157805, // group 7 band 3
        234.4047821555670, // group 7 band 4
        413.2123532117139, // group 7 band 5
        583.4942448986758, // group 7 band 6
        723.1626188903449, // group 7 band 7
        812.4107479403511  // group 7 band 8
    },
    {
        50.26168209703553, // group 8 band 1
        51.68607947024589, // group 8 band 2
        54.08606052191077, // group 8 band 3
        59.60911573317825, // group 8 band 4
        213.9038577022978, // group 8 band 5
        280.0828125915382, // group 8 band 6
        321.3354748837807, // group 8 band 7
        465.0872409987534  // group 8 band 8
    },
    {
        3.658884900516479, // group 9 band 1
        4.386512282709812, // group 9 band 2
        6.598633384651088, // group 9 band 3
        13.70067334645029, // group 9 band 4
        25.28131858801107, // group 9 band 5
        36.09000501594890, // group 9 band 6
        43.99485725811989, // group 9 band 7
        47.76971427517199  // group 9 band 8
    },
    {
        0.2311607986541810, // group 10 band 1
        0.2792642566974985, // group 10 band 2
        0.4285853026610446, // group 10 band 3
        0.9231299170458391, // group 10 band 4
        1.759261350905440,  // group 10 band 5
        2.568007049966189,  // group 10 band 6
        3.176657778709061,  // group 10 band 7
        3.473741688650271   // group 10 band 8
    }};
} // namespace benchmarkData

//---------------------------------------------------------------------------//

// declaration
bool checkData(rtt_dsxx::ScalarUnitTest &ut, SP_Goo spGandOpacity);

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);

  bool itPassed;

  // get the ipcress file name, and create the ipcress file
  string ipcressFileName =
      ut.getTestSourcePath() + std::string("odfregression10.ipcress");
  std::shared_ptr<IpcressFile const> file;
  try {
    file.reset(new IpcressFile(ipcressFileName));
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressFile object for "
            << "file \"" << ipcressFileName << "\":" << excpt.what();
    FAILMSG(message.str());
    cout << "Aborting tests.";
    return ut.numFails;
  }

  //load the Ipcress ODFMG Opacity
  std::shared_ptr<IpcressOdfmgOpacity const> spGandOpacity;

  try {
    spGandOpacity.reset(new IpcressOdfmgOpacity(
        file, benchmarkData::matID, benchmarkData::model,
        benchmarkData::reaction, benchmarkData::numBands));

    ostringstream message;
    message << "Successfully read Ipcress file \"" << ipcressFileName << "\"."
            << endl;
    PASSMSG(message.str());
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object with "
            << "data from \"" << ipcressFileName << "\":" << excpt.what();
    FAILMSG(message.str());
    cout << "Aborting tests.";
    return ut.numFails;
  }

  // check the data
  itPassed = checkData(ut, spGandOpacity);

  if (itPassed) {
    cout << "Read the correct data from " << ipcressFileName << "." << endl;
  } else {
    cout << "Read incorrect data from " << ipcressFileName << "." << endl;
  }

  // pack and then unpack it
  std::vector<char> packed;

  packed = spGandOpacity->pack();

  std::shared_ptr<IpcressOdfmgOpacity const> spUnpackedGandOpacity;

  try {
    spUnpackedGandOpacity.reset(new IpcressOdfmgOpacity(packed));
    PASSMSG("Successfully unpacked IpcressOdfmgOpacity.");
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to unpack "
            << "data from \"" << ipcressFileName << "\":" << excpt.what();
    FAILMSG(message.str());
    cout << "Aborting tests.";
    return ut.numFails;
  }

  // check the unpacked data
  itPassed = checkData(ut, spUnpackedGandOpacity);

  if (itPassed) {
    cout << "Read the correct data from unpacked IpcressOdfmgOpacity." << endl;
  } else {
    cout << "Read incorrect data from unpacked IpcressOdfmgOpacity." << endl;
  }

  //make sure it won't unpack as something else
  itPassed = false;
  try {
    std::shared_ptr<IpcressMultigroupOpacity> opacity(
        new IpcressMultigroupOpacity(packed));
  } catch (rtt_dsxx::assertion const &err) {
    itPassed = true;
    ostringstream message;
    message << "Good, we caught the following assertion, \n" << err.what();
    PASSMSG(message.str());
  }

  if (!itPassed) {
    ostringstream msg;
    msg << "Failed to catch an illegal packing asserion "
        << "(odfmg should not unpack as mg).";
    FAILMSG(msg.str());
  }

  return ut.numFails;
}

//---------------------------------------------------------------------------//
bool checkData(rtt_dsxx::ScalarUnitTest &ut, SP_Goo spGandOpacity) {
  Require(spGandOpacity);

  rtt_cdi::OpacityModelType omt(spGandOpacity->getOpacityModelType());
  if (omt == rtt_cdi::IPCRESS_TYPE)
    PASSMSG("OpacityModelType() returned expected value.");
  else
    FAILMSG("OpacityModelType() did not return the expected value.");

  std::string edp(spGandOpacity->getEnergyPolicyDescriptor());
  if (edp == std::string("odfmg"))
    PASSMSG("EDP = odfmg");
  else
    FAILMSG("EDP != odfmg");

  if (!spGandOpacity->data_in_tabular_form())
    ITFAILS;

  rtt_cdi::Model om(spGandOpacity->getModelType());
  FAIL_IF_NOT(om == rtt_cdi::ROSSELAND);

  rtt_cdi::Reaction rt(spGandOpacity->getReactionType());
  FAIL_IF_NOT(rt == rtt_cdi::ABSORPTION);

  string const expectedIpcressFileName =
      ut.getTestSourcePath() + std::string("odfregression10.ipcress");
  std::string dataFilename(spGandOpacity->getDataFilename());
  if (dataFilename != expectedIpcressFileName)
    ITFAILS;

  std::string ddesc(spGandOpacity->getDataDescriptor());
  if (ddesc != std::string("Multigroup Rosseland Absorption"))
    ITFAILS;

  double const temperature = benchmarkData::temp;
  double const density = benchmarkData::dens;

  size_t const numBands = spGandOpacity->getNumBands();
  size_t const numGroups = spGandOpacity->getNumGroups();

  bool hasNotFailed = true;

  // this message should never happen, as the numbands is an input
  // parameter
  if (numBands != benchmarkData::numBands) {
    FAILMSG("Number of bands does not match.");
    cout << "Aborting test.";
    hasNotFailed = false;
    return hasNotFailed;
  } else {
    PASSMSG("Number of bands matches.");
  }

  // test the number of groups
  if (numGroups != benchmarkData::numGroups) {
    FAILMSG("Number of groups does not match.");
    cout << "Aborting test.";
    hasNotFailed = false;
    return hasNotFailed;
  } else {
    PASSMSG("Number of groups matches.");
  }

  bool itFails = false;

  // test group boundaries
  vec_d const groupBoundaries = spGandOpacity->getGroupBoundaries();
  for (size_t group = 0; group < numGroups; group++) {
    if (!soft_equiv(groupBoundaries[group],
                    benchmarkData::groupBoundaries[group])) {
      itFails = true;
      cout << "Mismatch in group boundaries for group " << group + 1 << endl;
    }
  }

  if (itFails) {
    FAILMSG("Group boundaries did not match.");
    hasNotFailed = false;
  } else {
    PASSMSG("Group boundaries did match.");
  }

  // test band boundaries
  vec_d const bandBoundaries = spGandOpacity->getBandBoundaries();

  itFails = false;

  for (size_t band = 0; band < numBands; band++) {
    if (!soft_equiv(bandBoundaries[band],
                    benchmarkData::bandBoundaries[band])) {
      itFails = true;
      cout << "Mismatch in band boundaries for band " << band + 1 << endl;
    }
  }

  if (itFails) {
    FAILMSG("Band boundaries did not match.");
    hasNotFailed = false;
  } else {
    PASSMSG("Band boundaries did match.");
  }

  // test opacities
  vec2_d multiBandOpacities = spGandOpacity->getOpacity(temperature, density);

  itFails = false;

  for (size_t group = 0; group < numGroups; group++) {
    for (size_t band = 0; band < numBands; band++) {
      if (!soft_equiv(multiBandOpacities[group][band],
                      benchmarkData::opacities[group][band])) {
        cout << "Mismatch in opacity for group " << group + 1 << "band "
             << band + 1 << endl;
        itFails = true;
      }
    }
  }

  if (itFails) {
    FAILMSG("Opacities did not match.");
    hasNotFailed = false;
  } else {
    PASSMSG("Opacities did match.");
  }

  // check alternate accessors:

  std::vector<double> vtemp(1, temperature);
  std::vector<std::vector<std::vector<double>>> mbovt(
      spGandOpacity->getOpacity(vtemp, density));

  std::vector<double> vdens(1, density);
  std::vector<std::vector<std::vector<double>>> mbovd(
      spGandOpacity->getOpacity(temperature, vdens));

  for (size_t i = 0; i < mbovt.size(); ++i)
    for (size_t j = 0; j < mbovt[i].size(); ++j)
      for (size_t k = 0; k < mbovt[i][j].size(); ++k)
        if (!soft_equiv(mbovt[i][j][k], mbovd[i][j][k]))
          ITFAILS;

  // blah

  std::vector<double> opacs(numGroups * numBands, 0.0);
  spGandOpacity->getOpacity(temperature, vdens.begin(), vdens.end(),
                            opacs.begin());
  // for( size_t i=0; i<numGroups; ++i)
  //     std::cout << opacs[i] << std::endl;

  spGandOpacity->getOpacity(vtemp.begin(), vtemp.end(), density, opacs.begin());

  spGandOpacity->getOpacity(vtemp.begin(), vtemp.end(), vdens.begin(),
                            vdens.end(), opacs.begin());

  return hasNotFailed;
}

//---------------------------------------------------------------------------//
// end of tIpcressOdfmgOpacity.cc
//---------------------------------------------------------------------------//
