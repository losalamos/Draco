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
#include <vector>

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
int const matID = 10001;

int const numGroups = 33;
double const groupBoundaries[numGroups + 1] = {
    1.0e-05, 0.0178, 0.0316, 0.0562, 0.075, 0.1,   0.133, 0.178, 0.237,
    0.316,   0.422,  0.562,  0.75,   1,     1.33,  1.78,  2.37,  3.16,
    4.22,    5.62,   7.5,    10,     13.3,  17.8,  23.7,  31.6,  42.2,
    56.2,    75.0,   100.0,  133.0,  178.0, 237.0, 300.0};

int const numBands = 8;
double const bandBoundaries[numBands + 1] = {0,    0.03, 0.13, 0.33, 0.63,
                                             0.83, 0.93, 0.98, 1};

double const temp = 0.1;
double const dens = 0.1;

double const opacities[numGroups][numBands] = {
    {
        0.14502658180416503, // group 1 band 1
        0.14582254476324005, // group 1 band 2
        0.14764361301594983, // group 1 band 3
        0.15064040874685866, // group 1 band 4
        0.15404691099952028, // group 1 band 5
        0.1562486891705428,  // group 1 band 6
        0.15735841244601859, // group 1 band 7
        0.15789710855936398  // group 1 band 8
    },
    {
        0.13570034869421219, // group 2 band 1
        0.13627261152815381, // group 2 band 2
        0.13759768816724363, // group 2 band 3
        0.13983711856654915, // group 2 band 4
        0.14221323501194666, // group 2 band 5
        0.14367191680575131, // group 2 band 6
        0.1444022628017117,  // group 2 band 7
        0.14474338731853667  // group 2 band 8
    },
    {0.12115944711068351, 0.12203336286952249, 0.12407122675409814,
     0.12756077698956336, 0.13132696056061524, 0.13367147621391529,
     0.13485455301066862, 0.13540923815893749},
    {0.11137208851918494, 0.1119672809949552, 0.11335386870593078,
     0.11571685311721591, 0.11819751212018965, 0.11972376030717853,
     0.12049358180691831, 0.12085446650284505},
    {0.10010350556608978, 0.10078334261990415, 0.10237314729877321,
     0.10510140931213252, 0.10798998923843035, 0.10978018935123805,
     0.11068694821497969, 0.11111289842691802},
    {0.087610525460254593, 0.088349033918625239, 0.090084076964818904,
     0.093086817846837311, 0.096297700667840269, 0.098305351583954204,
     0.099327586518241195, 0.099809005160547415},
    {0.074027122373520077, 0.074807843322627451, 0.076653814431245942,
     0.079885948595909945, 0.08339089449874211, 0.085610097435852225,
     0.086748390573271369, 0.087286389032245323},
    {0.060648938972741415, 0.061396323369568719, 0.063175102016554505,
     0.066325894789003553, 0.069784272606174341, 0.072001583662114454,
     0.073147924239594334, 0.073691809696889182},
    {0.048068039418251336, 0.048748077295303199, 0.050379315824059037,
     0.053308629133221097, 0.056567878956158252, 0.058689755358438338,
     0.059797721228144064, 0.060325951152019337},
    {0.037047822588953873, 0.037624695172452406, 0.039019552528928879,
     0.041558275446532456, 0.044413191064043178, 0.046300938024635523,
     0.047297547766725526, 0.047775244252345833},
    {0.028135341070563077, 0.028592121227025539, 0.029703591653167185,
     0.031746001442515596, 0.034047457185530088, 0.035586727132836185,
     0.036407375231138232, 0.036802617574575321},
    {0.02114972649055584, 0.02150334310164739, 0.02236865834304155,
     0.023970177776068888, 0.025764220692679189, 0.026974341805525996,
     0.02762564306201707, 0.027940765897663294},
    {0.015870603325341428, 0.016138251741621003, 0.016795992068124172,
     0.018016892705174061, 0.019360814236208175, 0.020269278239807832,
     0.020761710672214599, 0.021000748681874345},
    {0.01193316608040048, 0.012134849231751503, 0.012632551503785212,
     0.013557252132390172, 0.014548182906230537, 0.015216815548603551,
     0.015581526805067821, 0.015759033179990414},
    {0.0089172073705403327, 0.0090729229360275668, 0.0094599905150222555,
     0.010181976219818709, 0.010929720410013404, 0.01143430907728344,
     0.011712532485101476, 0.011848575289198752},
    {0.0066970933301586785, 0.0068140632774088827, 0.0071065265390602246,
     0.0076517483836983828, 0.0081915183507141488, 0.0085534858387927501,
     0.0087549934489551656, 0.0088539034665366736},
    {0.0050230903753644843, 0.0051131470545129499, 0.005340221651820967,
     0.0057638211032714284, 0.0061613900231417265, 0.0064261670010764315,
     0.0065758403702597503, 0.0066497918950300904},
    {0.0037616270829508684, 0.0038313376607828327, 0.0040086061104911376,
     0.0043384880773063269, 0.0046301797609254484, 0.0048218935952957058,
     0.004932341293988605, 0.0049873909638595672},
    {0.0028246542601426576, 0.0028783338709440082, 0.003015506324464374,
     0.0032683972057042068, 0.0034796778576495924, 0.0036152515318590753,
     0.0036948384458916402, 0.0037349075202015142},
    {0.0021168501747625581, 0.0021591173106445759, 0.0022670844185386004,
     0.0024632216854649437, 0.0026199742790630007, 0.0027169962496540537,
     0.0027749842852847857, 0.0028045782564443968},
    {0.001587802326442275, 0.0016209834362758365, 0.0017047313174536169,
     0.00185390521745766, 0.001970121876744735, 0.002038964504467201,
     0.0020803290452794053, 0.0021017307878895431},
    {0.0011939712360112364, 0.0012200132647626029, 0.001284146476646946,
     0.001396330354382799, 0.001482791569353877, 0.0015320412667904474,
     0.0015612269615605148, 0.0015764829612102554},
    {0.00089238877110191998, 0.0009133445509005684, 0.00096314100573734402,
     0.0010494896641027054, 0.0011160808222539004, 0.001153105403082406,
     0.0011743078302641713, 0.0011853925207239789},
    {0.00067034306723774922, 0.00068647889814438137, 0.00072352199552166085,
     0.00078735176768426917, 0.00083641940125356025, 0.00086339658128275338,
     0.00087827798303598682, 0.00088589700091121703},
    {0.00050292338672275534, 0.00051547690078502183, 0.00054350093924001046,
     0.00059177246477012306, 0.00062888429988433895, 0.00064924073090126002,
     0.00066014923780426192, 0.00066548213638039505},
    {0.0003767431926875607, 0.00038638107270285533, 0.00040752760469780754,
     0.00044397857107297882, 0.0004720029066033414, 0.00048737390237675425,
     0.00049549759438501941, 0.0004992321019207577},
    {0.00028298325226629791, 0.00029023098585989669, 0.00030599094127397631,
     0.00033312410840196217, 0.00035393719108350825, 0.00036533768110107146,
     0.0003713360034017286, 0.0003739488556545859},
    {0.00021213624380454874, 0.00021762943543764482, 0.00022954082016345419,
     0.00025006606260452087, 0.0002658145072946141, 0.00027444155647564545,
     0.0002789801672077822, 0.00028088486812172855},
    {0.00015914936473101503, 0.00016326673642753491, 0.00017218560788739954,
     0.00018754317667871539, 0.00019931745939195383, 0.00020575951635425722,
     0.00020914854836554599, 0.00021054682353953936},
    {0.00011968424269434958, 0.0001227598608452887, 0.00012941603517340581,
     0.0001408577843495581, 0.0001496228775162454, 0.00015440896076496639,
     0.00015692519491278985, 0.00015795749889416867},
    {8.9457108840905263e-05, 9.1803330706097694e-05, 9.6885497905762266e-05,
     0.0001056402911375915, 0.00011237924941591048, 0.00011606292966135192,
     0.00011799986034622698, 0.00011879391849621992},
    {6.7189061039601443e-05, 6.8924894364595681e-05, 7.2676123553510393e-05,
     7.9120319779981169e-05, 8.4076464062845618e-05, 8.6780250460690591e-05,
     8.8198098851826882e-05, 8.8778792038151139e-05},
    {5.3021638825917566e-05, 5.4171402154189534e-05, 5.6623806800741785e-05,
     6.0738819810976598e-05, 6.3823612684247563e-05, 6.5477804279030348e-05,
     6.6336827072318425e-05, 6.6687017086281412e-05}};

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
    message << "Successfully read Ipcress file \"" << ipcressFileName
            << "\".\n";
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
    cout << "Read the correct data from " << ipcressFileName << ".\n";
  } else {
    cout << "Read incorrect data from " << ipcressFileName << ".\n";
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
    cout << "Read the correct data from unpacked IpcressOdfmgOpacity.\n";
  } else {
    cout << "Read incorrect data from unpacked IpcressOdfmgOpacity.\n";
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

  FAIL_IF_NOT(spGandOpacity->data_in_tabular_form());

  rtt_cdi::Model om(spGandOpacity->getModelType());
  FAIL_IF_NOT(om == rtt_cdi::ROSSELAND);

  rtt_cdi::Reaction rt(spGandOpacity->getReactionType());
  FAIL_IF_NOT(rt == rtt_cdi::ABSORPTION);

  string const expectedIpcressFileName =
      ut.getTestSourcePath() + std::string("odfregression10.ipcress");
  std::string dataFilename(spGandOpacity->getDataFilename());
  FAIL_IF_NOT(dataFilename == expectedIpcressFileName);

  std::string ddesc(spGandOpacity->getDataDescriptor());
  FAIL_IF_NOT(ddesc == std::string("Multigroup Rosseland Absorption"));

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

  spGandOpacity->getOpacity(vtemp.begin(), vtemp.end(), density, opacs.begin());

  spGandOpacity->getOpacity(vtemp.begin(), vtemp.end(), vdens.begin(),
                            vdens.end(), opacs.begin());

  return hasNotFailed;
}

//---------------------------------------------------------------------------//
// end of tIpcressOdfmgOpacity.cc
//---------------------------------------------------------------------------//
