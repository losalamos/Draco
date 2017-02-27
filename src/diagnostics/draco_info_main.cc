//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/draco_info_main.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: tstScalarUnitTest.cc 6864 2012-11-08 01:34:45Z kellyt $
//---------------------------------------------------------------------------//

#include "draco_info.hh"
#include "ds++/Assert.hh"
#include "ds++/XGetopt.hh"
#include <iostream>

int main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;
  try {
    bool version(false);
    bool brief(false);
    rtt_diagnostics::DracoInfo di;

    // Preparing to parse command line arguments.
    rtt_dsxx::XGetopt::csmap long_options;
    long_options['b'] = "brief";
    long_options['h'] = "help";
    long_options['v'] = "version";
    std::map<char, std::string> help_strings;
    help_strings['b'] = "print a brief message.";
    help_strings['v'] = "print version information and exit.";
    help_strings['h'] = "print this message.";
    rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

    int c(0);
    while ((c = program_options()) != -1) {
      switch (c) {
      case 'v': // --version
        version = true;
        break;

      case 'b': // --brief
        brief = true;
        break;

      case 'h':
        std::cout << program_options.display_help("draco_info") << std::endl;
        return 0;
        break;

      default:
        break;
      }
    }

    if (version)
      cout << di.versionReport();
    else if (brief)
      cout << di.briefReport();
    else
      cout << di.fullReport();
  } catch (rtt_dsxx::assertion &err) {
    std::string msg = err.what();
    std::cout << "ERROR: While running " << argv[0] << ", " << err.what()
              << std::endl;
    ;
    return 1;
  } catch (std::exception &err) {
    std::cout << "ERROR: While running " << argv[0] << ", " << err.what()
              << std::endl;
    ;
    return 1;
  } catch (...) {
    std::cout << "ERROR: While running " << argv[0] << ", "
              << "An unknown C++ exception was thrown" << std::endl;
    ;
    return 1;
  }

  return 0;
}

//---------------------------------------------------------------------------//
// end of draco_info_main.cc
//---------------------------------------------------------------------------//
