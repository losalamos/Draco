//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   compton2/cskrwx.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Basic reader/writer of csk Compton files
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/XGetopt.hh"
#include "units/PhysicalConstants.hh"
#include "units/UnitSystem.hh"
//#include <cmath>
//#include <iomanip>
#include <iostream>
#include <string>
//#include <memory>
//#include <sstream>
//#include <vector>

using std::cout;
using std::endl;
//using std::ios;
using std::string;

using rtt_dsxx::soft_equiv;

//----------------------------------------------------------------------------//
/*!
 * \brief Basic reader of the csk ASCII file format
 *
 * Modification of this executable to make it more useful is encouraged.
 */
//----------------------------------------------------------------------------//
void read_csk_files(std::string const &basename) {
  // csk data base filename (csk ASCII format required)

  std::array<std::string, 2> inouts = {"in", "out"};
  std::array<std::string, 2> lins = {"lin", "nonlin"};

  // Normalization / unit change
  const double mtocm = 100.0;
  const double classical_electron_radius =
      mtocm * rtt_units::classicalElectronRadiusSI; // cm

  // Normalization constant for raw CSK data:
  // CSK to cross section: 2 * pi * classicalElectronRadius^2 / 4
  // cross section to opacity: Zbar_over_A * avogadrosNumber
  const double csk_opac_norm = 0.25 * 2 * rtt_units::PI *
                               classical_electron_radius *
                               classical_electron_radius * rtt_units::AVOGADRO;
  // this conversion is used to go from cm^2/mole to cm^2 (from opacity to micro xs)
  const double csk_xs_conv = 1. / rtt_units::AVOGADRO;

  // nonlinear scale is ~ 1.545956458889162e+26, but temperature-dependent

  Ensure(soft_equiv(csk_opac_norm, 2 * 0.037558, 1e-4));
  Ensure(soft_equiv(2 * 0.2003102 * csk_xs_conv, 6.6524587e-25, 1e-4));

  for (std::string lin : lins) {
    for (std::string inout : inouts) {

      std::string const filename = basename + '_' + inout + '_' + lin;
      cout << "Reading file: " << filename << endl;
    }
  }
}

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  // Process known command line arguments:
  rtt_dsxx::XGetopt::csmap long_options;
  long_options['h'] = "help";
  long_options['v'] = "version";
  std::map<char, std::string> help_strings;
  help_strings['h'] = "print this message.";
  help_strings['v'] = "print version information and exit.";
  rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

  std::string const helpstring(
      "\nUsage: cskrw [-hv] "
      "<csk_base_filename>\n¡Under active development!\n");

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'v': // --version
      cout << argv[0] << ": version " << rtt_dsxx::release() << endl;
      return 0;

    case 'h': // --help
      cout << argv[0] << ": version " << rtt_dsxx::release() << helpstring
           << endl;
      return 0;
    }
  }

  // Assume last command line argument is the name of the ipcress file.
  std::string const filename = string((argc > 1) ? argv[argc - 1] : "csk");

  try {
    read_csk_files(filename);
  } catch (rtt_dsxx::assertion &excpt) {
    cout << "While attempting to read csk file, " << excpt.what() << endl;
    return 1;
  }

  return 0;
}

//----------------------------------------------------------------------------//
// end of cskrw.cc
//----------------------------------------------------------------------------//
