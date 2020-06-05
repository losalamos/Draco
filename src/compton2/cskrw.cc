//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   compton2/cskrw.cc
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
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
//#include <memory>
//#include <sstream>
//#include <vector>

using std::cout;
using std::endl;
//using std::ios;
using std::string;

using rtt_dsxx::soft_equiv;

using UINT = uint64_t;
using FP = double;

//----------------------------------------------------------------------------//
/*!
 * \brief Basic reader of the csk ASCII file format
 *
 * Modification of this executable to make it more useful is encouraged.
 *
 * What should this do? Read ASCII and do ... ?
 * Different than Converter++, which converts ASCII to binary?
 * Different than funcitonality to interpolate in T and print out a few files?
 *
 * Will be useful in testing ASCII reader and renormalization!!
 */
//----------------------------------------------------------------------------//
void read_csk_files(std::string const &basename, int verbosity) {
    // csk data base filename (csk ASCII format required)

#if 1
    std::array<std::string, 2> inouts = {"in", "out"};
    std::array<std::string, 2> lins = {"lin", "nonlin"};
#else
    std::array<std::string, 1> inouts = {"in"};
    std::array<std::string, 1> lins = {"lin"};
#endif

    // Normalization / unit change
    const FP mtocm = 100.0;
    const FP classical_electron_radius =
        mtocm * rtt_units::classicalElectronRadiusSI; // cm

    // Normalization constant for raw CSK data:
    // CSK to cross section: 2 * pi * classicalElectronRadius^2 / 4
    // cross section to opacity: Zbar_over_A * avogadrosNumber
    const FP csk_opac_norm = 0.25 * 2 * rtt_units::PI *
        classical_electron_radius *
        classical_electron_radius * rtt_units::AVOGADRO;
    // this conversion is used to go from cm^2/mole to cm^2 (from opacity to micro xs)
    const FP csk_xs_conv = 1. / rtt_units::AVOGADRO;

    // Normalization for T and E (keV; value used by CSK)
    const FP mec2 = 510.9989461;

    // nonlinear scale is ~ 1.545956458889162e+26, but temperature-dependent

    Ensure(soft_equiv(csk_opac_norm, 2 * 0.037558, 1e-4));
    Ensure(soft_equiv(2 * 0.2003102 * csk_xs_conv, 6.6524587e-25, 1e-4));

    const UINT numEvals = lins.size() * inouts.size();
    UINT counter = 0;
    std::vector<FP> data;

    for (std::string lin : lins) {
        for (std::string inout : inouts) {

            ++counter;
            UINT eval = counter - 1U;

            std::string const filename = basename + '_' + inout + '_' + lin;
            cout << "Reading file: " << filename << endl;

            std::ifstream f(filename);
            Insist(f.is_open(), "Unable to open " + filename);

            // set effective renormalization constant
            // TODO: Check if this is right; change for nonlinear
            const FP renorm = csk_opac_norm; // cm^2/mole
            if (verbosity > 0)
                std::cout << "renorm: " << renorm << '\n';

            // Should later encapsulate read and store
            {
                // Line 1: sizes
                UINT numTbreakpoints = 0;
                UINT numTs = 0;
                UINT numGroups = 0;
                UINT numLegMoments = 0;
                f >> numTbreakpoints >> numTs >> numGroups >> numLegMoments;
                if (verbosity > 0)
                    std::cout << "CSK sizes: " << numTbreakpoints << ' ' << numTs << ' ' << numGroups << ' ' << numLegMoments << '\n';

                // Line 2: T breakpoints
                std::vector<FP> Tbreakpoints(numTbreakpoints, 0.0);
                if (verbosity > 0)
                    std::cout << "Temperature breakpoints (keV):";
                for (UINT i=0; i < numTbreakpoints; ++i)
                {
                    f >> Tbreakpoints[i];
                    // scale to keV
                    Tbreakpoints[i] *= mec2;
                    if (verbosity > 0)
                        std::cout << ' ' << std::setprecision(4) << Tbreakpoints[i];
                }
                if (verbosity > 0)
                    std::cout << '\n';

                // Line 3: Group bounds
                std::vector<FP> groupBdrs(numGroups + 1U, 0.0);
                if (verbosity > 0)
                    std::cout << "Group boundaries (keV):";
                for (UINT g=0; g < numGroups+1U; ++g)
                {
                    f >> groupBdrs[g];
                    // scale to keV
                    groupBdrs[g] *= mec2;
                    if (verbosity > 0)
                        std::cout << ' ' << std::setprecision(4) << groupBdrs[g];
                }
                if (verbosity > 0)
                    std::cout << '\n';

                // Remaining lines: Temperatures and MG data
                std::vector<FP> Ts(numTs, 0.0);
                // data is 1D array indexed [eval, legmoment, groupfrom, groupto]
                if (!data.size())
                {
                    data.resize(numEvals * numLegMoments * numGroups * numGroups);
                }
                // Format:
                // T
                // gto gfrom moment0 [moment1 moment2 ...]
                // <blankline>
                // tmp storage
                std::string line;
                for (UINT iT=0; iT < numTs; ++iT)
                {
                    f >> Ts[iT];
                    // scale to keV
                    Ts[iT] *= mec2;
                    if (verbosity > 0)
                        std::cout << "T(keV) " << Ts[iT] << '\n';

                    UINT gfrom = 0U;
                    UINT gto = 0U;
                    bool finished = false;
                    while (!finished)
                    {
                        // Read one line
                        f >> gfrom >> gto;
                        if (verbosity > 1)
                            std::cout << "  " << gfrom << ' ' << gto;
                        // 1-based to 0-based
                        gfrom -= 1U;
                        gto -= 1U;

                        // Read xs
                        for (UINT iL=0; iL < numLegMoments; ++iL)
                        {
                            // Read value
                            FP val;
                            f >> val;
                            val *= renorm;

                            // Put in 1D data vector
                            const UINT loc = gto + numGroups * (gfrom + numGroups * (iL + numLegMoments * eval));
                            Require(loc < data.size());
                            data[loc] = val;
                            //std::cout << ' ' << std::setprecision(2) << val / renorm;
                            if (verbosity > 1)
                                std::cout << ' ' << std::setprecision(2) << val;
                        }
                        if (verbosity > 1)
                            std::cout << '\n';

                        // Look-ahead to next line to see if a blank line,
                        // signalling end of temperature data
                        // TODO: Make this not use getline
                        std::getline(f, line);
                        int next = f.get();
                        if (next == '\n' || next == EOF)
                        {
                            finished = true;
                            if (verbosity > 1)
                                std::cout << '\n';
                        }
                        else
                        {
                            f.putback(char(next));
                        }
                    }
                }
            }

            f.close();
            if (verbosity > 0)
                std::cout << '\n';
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
        int verbosity = 0;
        verbosity = 1;
        //verbosity = 2;
        read_csk_files(filename, verbosity);
    } catch (rtt_dsxx::assertion &excpt) {
        cout << "While attempting to read csk file, " << excpt.what() << endl;
        return 1;
    }

    return 0;
}

//----------------------------------------------------------------------------//
// end of cskrw.cc
//----------------------------------------------------------------------------//
