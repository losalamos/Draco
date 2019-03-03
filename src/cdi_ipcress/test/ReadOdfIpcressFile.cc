//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/test/ReadOdfIpcressFile.cc
 * \author Seth R. Johnson
 * \date   Thu July 10 2008
 * \note   Copyright (C) 2008-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "cdi_ipcress_test.hh"
#include "cdi_ipcress/IpcressOdfmgOpacity.hh"
#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"
#include <cstdlib> // atof

using rtt_cdi_ipcress::IpcressFile;
using rtt_cdi_ipcress::IpcressOdfmgOpacity;

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::istringstream;
using std::ostringstream;
using std::string;

typedef std::shared_ptr<const IpcressOdfmgOpacity> SP_Goo;
typedef std::vector<double> vec_d;
typedef std::vector<std::vector<double>> vec2_d;

void printGrid(SP_Goo spGandOpacity);
void askTempDens(double &temperature, double &density, bool is_unittest);
void analyzeData(SP_Goo spGandOpacity);
void collapseOpacities(SP_Goo spGandOpacity, double temperature,
                       double density);
void printCData(SP_Goo spGandOpacity, double temperature, double density);
void printData(SP_Goo spGandOpacity, double temperature, double density);
void printTable(SP_Goo spGandOpacity, double temperature, double density);

//---------------------------------------------------------------------------//
void printSyntax() {
  cerr << "syntax: ReadOdfIpcressFile [--help] [--bands n | --mg] "
       << "[--analyze | --printc | --printtable | --collapse] [--model "
       << "(r[osseland] | p[lanck])] [--reaction (t[otal] | a[bsorption]"
       << " | s[cattering])] [-d density -t temperature] fileName" << endl;

  // Undocumented:
  // ReadOdfIpcressFile --unittest
}
//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  bool is_unittest(false);

  // Known command line arguments:
  // rtt_dsxx::XGetopt::csmap long_options;
  // long_options['h'] = "help";
  // long_options['v'] = "version";
  // long_options['u'] = "unittest";
  // long_options['m'] = "model:";
  // long_options['r'] = "reaction:";
  // long_options['g'] = "mg";
  // long_options['b'] = "bands";
  // long_options['a'] = "analyze";
  // long_options['p'] = "printc";
  // long_options['c'] = "collapse";
  // long_options['i'] = "printtable";
  // long_options['c'] = "density:";
  // long_options['t'] = "temperature:";
  // std::map<char,std::string> help_strings;
  // help_strings['h'] = "print this message.";
  // help_strings['v'] = "print version information and exit.";
  // help_strings['u'] = "unittest mode";
  // help_strings['m'] = "model value = r[osseland] | p[lanck]";
  // help_strings['r'] = "reaction value = t[otal] | a[bsorption] | s[cattering]";
  // help_strings['g'] = "This option cannot be combined with --bands";
  // help_strings['b'] = "This options cannot be combined with --mg";
  // help_strings['a'] = "This option cannot be combined with printc, printtable or collapse";
  // help_strings['p'] = "This option cannot be combined with analyze, printtable or collapse";
  // help_strings['c'] = "This option cannot be combined with analyze, printc or collapse";
  // help_strings['i'] = "This option cannot be combined with analyze, printc, printtable or collapse";
  // help_strings['c'] = "density";
  // help_strings['t'] = "temperature";
  // rtt_dsxx::XGetopt program_options( argc, argv, long_options, help_strings );

  if (argc <= 1) {
    cerr << "Must have at least one input argument (the ipcress file)." << endl;
    printSyntax();
    return ut.numFails++;
  } else if (argc == 2) {
    rtt_dsxx::XGetopt program_options(argc, argv, "hu");
    int c(0);
    while ((c = program_options()) != -1) {
      switch (c) {
      case 'h': // --help
        printSyntax();
        return 0;
        break;

      case 'u': // --unittest
        is_unittest = true;
        break;

      default:
        return 0; // nothing to do.
      }
    }
  }

  // get the ipcress file name, and create the ipcress file
  string ipcressFileName;
  if (!is_unittest)
    ipcressFileName = argv[argc - 1];
  else
    ipcressFileName =
        ut.getTestSourcePath() + std::string("odfregression10.ipcress");

  std::shared_ptr<const IpcressFile> file;
  try {
    file.reset(new IpcressFile(ipcressFileName));
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create SP to new IpcressFile object for "
            << "file \"" << ipcressFileName << "\": " << excpt.what();
    FAILMSG(message.str());
    cerr << "Aborting program.";
    return 1;
  }

  // data needed to load the Ipcress file properly
  int numBands = 0;
  rtt_cdi::Model model = rtt_cdi::ROSSELAND;
  rtt_cdi::Reaction reaction = rtt_cdi::ABSORPTION;
  int matID = 19000;
  int actionToTake = 0; // 0 for user input, 1 for analyze, 2 for print c

  double temperature = 0;
  double density = 0;

  // loop on all arguments except the first (program name) and last (input file
  // name)
  for (int arg = 1; arg < argc - 1; arg++) {
    rtt_dsxx::XGetopt program_options(argc, argv, "hdtmrgbapci");

    int c(0);
    while ((c = program_options()) != -1) {
      switch (c) {
      case 'h':
        printSyntax();
        return 0;
        break;

      case 'd':
        // density = std::stod(program_options.get_option_value()); // C++11
        density = atof(program_options.get_option_value().c_str());
        cerr << "Using density of " << density << endl;
        break;

      case 't':
        temperature = atof(program_options.get_option_value().c_str());
        cerr << "Using temperature of " << temperature << endl;
        break;

      case 'm': // --model
      {
        std::string const parsed_model = program_options.get_option_value();
        if (parsed_model[0] == 'r')
          model = rtt_cdi::ROSSELAND;
        else if (parsed_model[0] == 'p')
          model = rtt_cdi::PLANCK;
        else
          cerr << "Unrecognized model option '" << parsed_model << "'\n"
               << "Defaulting to rosseland" << endl;
        break;
      }

      case 'r': // --reaction
      {
        std::string const parsed_reaction = program_options.get_option_value();
        if (parsed_reaction[0] == 'a')
          reaction = rtt_cdi::ABSORPTION;
        else if (parsed_reaction[0] == 's')
          reaction = rtt_cdi::SCATTERING;
        else if (parsed_reaction[0] == 't')
          reaction = rtt_cdi::TOTAL;
        else
          cerr << "Unrecognized model option '" << parsed_reaction << "'\n"
               << "Defaulting to rosseland" << endl;
        break;
      }
      case 'g': // --mg
        numBands = 1;
        cerr << "Using " << numBands << " bands (multigroup file)" << endl;
        break;

      case 'b': // --bands
        numBands = atoi(program_options.get_option_value().c_str());
        cerr << "Using " << numBands << " bands" << endl;
        break;

      case 'a': // --analyze
        actionToTake = 1;
        break;

      case 'p': // --printc
        actionToTake = 2;
        break;

      case 'c': // --collapse
        actionToTake = 3;
        break;

      case 'i': // --printtable
        actionToTake = 4;
        break;

      default:
        std::vector<std::string> const &ua =
            program_options.get_unmatched_arguments();
        cerr << "Error: " << ua.size() << " unrecognized options:\n";
        for (size_t index = 0; index < ua.size(); ++index)
          std::cout << "\n\t" << ua[index];
        return 0;
        break;
      }
    }
  }

  //print the model that we're using
  if (model == rtt_cdi::ROSSELAND)
    cerr << "Using ROSSELAND weighting" << endl;
  else
    cerr << "Using PLANCK weighting" << endl;

  //print the cross section that we're using
  if (reaction == rtt_cdi::TOTAL)
    cerr << "Using TOTAL reaction" << endl;
  else if (reaction == rtt_cdi::ABSORPTION)
    cerr << "Using ABSORPTION reaction" << endl;
  else
    cerr << "Using SCATTERING reaction" << endl;

  //ask the user for the number of bands
  if (is_unittest)
    numBands = 1;
  while (numBands == 0) {
    cout << "Enter the number of bands (use 1 for multigroup file): ";
    cin >> numBands;
  }
  Insist(numBands > 0, "Must have a positive number of bands.");

  //load the Ipcress ODFMG Opacity
  std::shared_ptr<const IpcressOdfmgOpacity> spGandOpacity;

  spGandOpacity.reset(
      new IpcressOdfmgOpacity(file, matID, model, reaction, numBands));

  cerr << "Successfully read Ipcress file \"" << ipcressFileName << "\"."
       << endl;

  switch (actionToTake) {
  case 0:
    if (rtt_dsxx::soft_equiv(temperature, 0.0) ||
        rtt_dsxx::soft_equiv(density, 0.0)) {
      printGrid(spGandOpacity);
      askTempDens(temperature, density, is_unittest);
    }
    printData(spGandOpacity, temperature, density);
    break;
  case 1:
    analyzeData(spGandOpacity);
    break;
  case 2:
    if (rtt_dsxx::soft_equiv(temperature, 0.0) ||
        rtt_dsxx::soft_equiv(density, 0.0)) {
      printGrid(spGandOpacity);
      askTempDens(temperature, density, is_unittest);
    }
    printCData(spGandOpacity, temperature, density);
    break;
  case 3:
    if (rtt_dsxx::soft_equiv(temperature, 0.0) ||
        rtt_dsxx::soft_equiv(density, 0.0)) {
      printGrid(spGandOpacity);
      askTempDens(temperature, density, is_unittest);
    }
    collapseOpacities(spGandOpacity, temperature, density);
    break;
  case 4:
    if (rtt_dsxx::soft_equiv(temperature, 0.0) ||
        rtt_dsxx::soft_equiv(density, 0.0)) {
      printGrid(spGandOpacity);
      askTempDens(temperature, density, is_unittest);
    }
    printTable(spGandOpacity, temperature, density);

    break;
  }

  // Check status of spGandolfOpacity if this is a unit test
  if (is_unittest) {
    std::cout << "\nChecking a few values for the unit test...\n" << std::endl;
    std::vector<std::vector<double>> opac =
        spGandOpacity->getOpacity(5.0, 0.05);
    if (opac.size() != 80)
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(opac[0][0], 2128.526464249052))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(opac[10][0], 221.5324065688365))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(opac[20][0], 12.04514705449304))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(opac[30][0], 0.9751198562573208))
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(opac[40][0], 0.3344851514293186))
      ITFAILS;
    if (ut.numFails != 0) {
      std::cout.precision(16);
      std::cout << "opac[0][0]  = " << opac[0][0] << std::endl;
      std::cout << "opac[10][0] = " << opac[10][0] << std::endl;
      std::cout << "opac[20][0] = " << opac[20][0] << std::endl;
      std::cout << "opac[30][0] = " << opac[30][0] << std::endl;
      std::cout << "opac[40][0] = " << opac[40][0] << std::endl;
    }
    std::vector<double> grp_bnds = spGandOpacity->getGroupBoundaries();
    if (grp_bnds.size() != 81)
      ITFAILS;
    if (!rtt_dsxx::soft_equiv(grp_bnds[0], 0.01) ||
        !rtt_dsxx::soft_equiv(grp_bnds[80], 100.0))
      ITFAILS;
    if (spGandOpacity->getNumTemperatures() != 6)
      ITFAILS;
    if (spGandOpacity->getNumDensities() != 4)
      ITFAILS;
    if (ut.numFails == 0)
      PASSMSG("Successfully extracted data from odfregression10.ipcress.");
    else
      FAILMSG("Failed to extract expected data from odfregression10.ipcress.");
  }

  cerr << "\nFinished." << endl;
  return 0;
}
//---------------------------------------------------------------------------//
void printGrid(SP_Goo spGandOpacity) {
  const vec_d temperatures = spGandOpacity->getTemperatureGrid();
  const vec_d densities = spGandOpacity->getDensityGrid();

  cout << "Temperature grid: ";
  for (size_t i = 0; i < temperatures.size(); i++)
    cout << temperatures[i] << " ";
  cout << endl;

  cout << "Density grid: ";
  for (size_t i = 0; i < densities.size(); i++)
    cout << densities[i] << " ";
  cout << endl;
}
//---------------------------------------------------------------------------//
void askTempDens(double &temperature, double &density, bool is_unittest) {
  if (is_unittest) {
    temperature = 5.0;
    density = 0.05;
    return;
  }
  while (temperature <= 0 || density <= 0) {
    cout << "Enter the temperature to analyze: ";
    cin >> temperature;

    cout << "Enter the density to analyze    : ";
    cin >> density;
  }
}
//---------------------------------------------------------------------------//
void analyzeData(SP_Goo spGandOpacity) {
  Require(spGandOpacity);

  const size_t numBands = spGandOpacity->getNumBands();
  const size_t numGroups = spGandOpacity->getNumGroups();

  const vec_d temperatures = spGandOpacity->getTemperatureGrid();
  const vec_d densities = spGandOpacity->getDensityGrid();

  cout << "Temperature\tDensity\tGroup opacity range\n";

  //loop over all stored temperatures and densities
  for (size_t t = 0; t < temperatures.size(); t++) {
    for (size_t d = 0; d < densities.size(); d++) {
      vec2_d multiBandOpacities =
          spGandOpacity->getOpacity(temperatures[t], densities[d]);
      cout << temperatures[t] << "\t" << densities[d] << "\t";

      for (size_t group = 0; group < numGroups; group++) {
        cout << multiBandOpacities[group][numBands - 1] /
                    multiBandOpacities[group][0]
             << "\t";
      }
      cout << endl;
    }
  }
}
//---------------------------------------------------------------------------//
void collapseOpacities(SP_Goo spGandOpacity, double temperature,
                       double density) {
  Require(spGandOpacity);

  const size_t numBands = spGandOpacity->getNumBands();
  const size_t numGroups = spGandOpacity->getNumGroups();
  const rtt_cdi::Model model = spGandOpacity->getModelType();

  cout << "=============================================" << endl;
  cout << "Printing collapsed opacities at " << temperature << " keV, "
       << "rho = " << density << endl;
  cout << "=============================================" << endl;

  const vec_d groupBoundaries = spGandOpacity->getGroupBoundaries();
  const vec_d bandBoundaries = spGandOpacity->getBandBoundaries();

  vec_d bandWidths(numBands, 0.0);

  //calculate band widths
  Check(static_cast<size_t>(numBands + 1) == bandBoundaries.size());

  for (size_t band = 0; band < numBands; band++) {
    bandWidths[band] = bandBoundaries[band + 1] - bandBoundaries[band];
  }

  //loop over groups
  vec2_d multiBandOpacities = spGandOpacity->getOpacity(temperature, density);

  for (size_t group = 0; group < numGroups; group++) {
    double collapsedOpacity = 0;

    // harmonic average for rosseland
    if (model == rtt_cdi::ROSSELAND) {
      for (size_t band = numBands - 1; band != 0; band--) {
        collapsedOpacity += bandWidths[band] / multiBandOpacities[group][band];
      }
      if (!rtt_dsxx::soft_equiv(collapsedOpacity, 0.0))
        collapsedOpacity = 1 / collapsedOpacity;
    } else // arithmetic average for planckian
    {
      for (size_t band = 0; band < numBands; band++) {
        collapsedOpacity += bandWidths[band] * multiBandOpacities[group][band];
      }
    }

    printf("%4zu\t%.6f", group + 1, collapsedOpacity);
    cout << endl;
  }
}

//---------------------------------------------------------------------------//
void printTable(SP_Goo spGandOpacity, double temperature, double density) {
  Require(spGandOpacity);

  const size_t numBands = spGandOpacity->getNumBands();
  const size_t numGroups = spGandOpacity->getNumGroups();

  cout << "Temperature:\t" << temperature << "\tDensity:\t" << density << endl;

  // prsize_t group boundaries
  cout << "numGroups:\t" << numGroups << "\tnumBands:\t" << numBands << endl;

  const vec_d groupBoundaries = spGandOpacity->getGroupBoundaries();

  for (size_t i = 0; i < groupBoundaries.size(); i++) {
    printf("%.6g\t", groupBoundaries[i]);
  }
  cout << endl;

  const vec_d bandBoundaries = spGandOpacity->getBandBoundaries();

  for (size_t i = 0; i < bandBoundaries.size(); i++) {
    printf("%.6g\t", bandBoundaries[i]);
  }
  cout << endl;

  // print opacity data
  vec2_d multiBandOpacities = spGandOpacity->getOpacity(temperature, density);

  for (size_t group = 0; group < numGroups; group++) {
    // print data for each band
    for (size_t band = 0; band < numBands; band++) {
      printf("%.16g\t", multiBandOpacities[group][band]);
    }

    cout << endl;
  }
  cout << endl;
}

//---------------------------------------------------------------------------//
void printCData(SP_Goo spGandOpacity, double temperature, double density) {
  Require(spGandOpacity);

  const size_t numBands = spGandOpacity->getNumBands();
  const size_t numGroups = spGandOpacity->getNumGroups();

  cout << "=============================================" << endl;
  cout << "Printing data at " << temperature << " keV, "
       << "rho = " << density << endl;
  cout << "=============================================" << endl;

  // prsize_t group boundaries
  cout << "const int numGroups = " << numGroups << ";" << endl;
  cout << "const double groupBoundaries[numGroups + 1] = {" << endl;

  const vec_d groupBoundaries = spGandOpacity->getGroupBoundaries();

  for (size_t i = 0; i < groupBoundaries.size(); i++) {
    printf("\t%.6g", groupBoundaries[i]);
    if (i != groupBoundaries.size() - 1)
      cout << ",";

    cout << endl;
  }
  cout << "};" << endl;

  // print band boundaries

  cout << "const int numBands = " << numBands << ";" << endl;
  cout << "const double bandBoundaries[numBands + 1] = {" << endl;

  const vec_d bandBoundaries = spGandOpacity->getBandBoundaries();

  for (size_t i = 0; i < bandBoundaries.size(); i++) {
    printf("\t%.6g", bandBoundaries[i]);
    if (i != bandBoundaries.size() - 1)
      cout << ",";

    cout << endl;
  }
  cout << "};" << endl;

  // print opacity data
  cout << "const double opacities[numGroups][numBands] = {\n";

  vec2_d multiBandOpacities = spGandOpacity->getOpacity(temperature, density);

  for (size_t group = 0; group < numGroups; group++) {
    cout << "{\n";
    // print data for each band
    for (size_t band = 0; band < numBands; band++) {
      printf("\t%#25.16g", multiBandOpacities[group][band]);

      if (band != numBands - 1)
        cout << ",";

      printf("\t\t// group %zu band %zu", group + 1, band + 1);

      cout << "\n";
    }

    cout << "}";
    if (group != numGroups - 1)
      cout << ",";

    cout << "\n";
  }
  cout << "};" << endl;
}

//---------------------------------------------------------------------------//
void printData(SP_Goo spGandOpacity, double temperature, double density) {
  Require(spGandOpacity);

  const size_t numBands = spGandOpacity->getNumBands();
  const size_t numGroups = spGandOpacity->getNumGroups();

  const vec_d groupBoundaries = spGandOpacity->getGroupBoundaries();
  const vec_d bandBoundaries = spGandOpacity->getBandBoundaries();

  cout << "=============================================\n";
  cout << "Printing band data at " << temperature << " keV, "
       << "rho = " << density << "\n";

  vec2_d multiBandOpacities = spGandOpacity->getOpacity(temperature, density);

  size_t maxGroup = 0;
  double maxRatio = 0.0;

  for (size_t group = 0; group < numGroups; group++) {
    double currentRatio = 0.0;

    cout << "=== Group " << group + 1 << " has energy range ["
         << groupBoundaries[group] << "," << groupBoundaries[group + 1] << "]\n"
         << "Group Band  Width        Opacity  Ratio to first\n";

    // print data for each band
    for (size_t band = 0; band < numBands; band++) {
      currentRatio =
          multiBandOpacities[group][band] / multiBandOpacities[group][0];

      printf("%5zu %4zu %6.3f %#14.6G  %14.3f\n", group + 1, band + 1,
             bandBoundaries[band + 1] - bandBoundaries[band],
             multiBandOpacities[group][band], currentRatio);
    }

    // compare result for last band
    if (currentRatio > maxRatio) {
      maxRatio = currentRatio;
      maxGroup = group;
    }
  }

  cout << "=============================================\n"
       << "At " << temperature << " keV, "
       << "rho = " << density << "\n"
       << "Best odf was in group " << maxGroup << " which had a high-to-low "
       << "ratio of " << maxRatio << ".\n"
       << "=============================================" << endl;
}

//---------------------------------------------------------------------------//
// end of ReadOdfIpcressFile.cc
//---------------------------------------------------------------------------//
