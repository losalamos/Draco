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
#include "cdi/CDI.hh" //planck integral
#include <cmath>
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

struct Dense_Compton_Data
{
    UINT numEvals;
    UINT numTs;
    UINT numGroups;
    UINT numLegMoments;
    std::vector<FP> groupBdrs;
    std::vector<FP> Ts;
    // [eval, moment, T, gfrom, gto]
    std::vector<FP> data;
    // [eval, moment, T, gfrom, gto]
    std::vector<FP> derivs;

    void resize(UINT numfiles, std::string filename);
    void read_from_file(UINT eval, std::string filename, bool isnonlin);
    void print_contents(int verbosity, int precision);
};

// resize data and set sizes
// numfiles is number of Compton files; filename is one such file
void Dense_Compton_Data::resize(UINT numfiles, std::string filename)
{
    // Reserve space for nl difference
    UINT numderived = (numfiles >= 4) ? 1 : 0;
    numEvals = numfiles + numderived;

    // Read first line of filename to get sizes
    std::ifstream f(filename);
    Insist(f.is_open(), "Unable to open " + filename);

    // Line 1: sizes
    UINT numTbreakpoints = 0; // not used
    numTs = 0;
    numGroups = 0;
    numLegMoments = 0;
    f >> numTbreakpoints >> numTs >> numGroups >> numLegMoments;

    // Set vector lengths
    groupBdrs.resize(numGroups+1, 0.0);
    Ts.resize(numTs, 0.0);
    UINT sz = numEvals * numLegMoments * numTs * numGroups * numGroups;
    data.resize(sz, 0.0);
    derivs.resize(sz, 0.0);

    f.close();
}

// read the entire contents of one file
void Dense_Compton_Data::read_from_file(UINT eval, std::string filename, bool isnonlin)
{
    Insist(eval < numEvals, "eval must be < numEvals");
    std::ifstream f(filename);
    Insist(f.is_open(), "Unable to open " + filename);

    // Normalization / unit change
    const FP mtocm = 100.0;
    const FP classical_electron_radius =
        mtocm * rtt_units::classicalElectronRadiusSI; // cm

    // Normalization constants for raw CSK data:
    // CSK to cross section: 2 * pi * classicalElectronRadius^2 / 4
    // cross section to opacity: Zbar_over_A * avogadrosNumber
    // opacity to micro xs: 1/avogadrosNumber
    //
    // convert from CSK data to opacity (cm^2/mole)
    const FP csk_opac_norm = 0.25 * 2 * rtt_units::PI *
        classical_electron_radius *
        classical_electron_radius * rtt_units::AVOGADRO;
    Ensure(soft_equiv(csk_opac_norm, 2 * 0.037558, 1e-4));
    //
    // convert from opacity (cm^2/mole) to micro xs (cm^2)
    //const FP csk_xs_conv = 1. / rtt_units::AVOGADRO;
    //Ensure(soft_equiv(2 * 0.2003102 * csk_xs_conv, 6.6524587e-25, 1e-4));

    // VALUES USED IN CSK'S PHYSICAL_CONSTANTS.HH
    // Normalization for energy and nonlinear terms
    const FP mec2 = 510.998;
    const FP hplanck = 8.6173303e-8;
    const FP cspeed = 299.79245800;

    // set effective renorm base multipliers
    const FP basescale = csk_opac_norm; // cm^2/mole
    const FP nlbase = (4.0/9.0) * 2.0 / (hplanck * hplanck * hplanck * cspeed * cspeed) * (mec2 * mec2 * mec2);

    // Line 1: sizes
    UINT numTbreakpoints = 0;
    UINT numTs_check = 0;
    UINT numGroups_check = 0;
    UINT numLegMoments_check = 0;
    f >> numTbreakpoints >> numTs_check >> numGroups_check >> numLegMoments_check;
    Check(numTs == numTs_check);
    Check(numGroups == numGroups_check);
    Check(numLegMoments == numLegMoments_check);

    // Line 2: T breakpoints (unused)
    FP throwaway;
    for (UINT i=0; i < numTbreakpoints; ++i)
    {
        f >> throwaway;
    }

    // Line 3: Group bounds
    for (UINT g=0; g < numGroups+1U; ++g)
    {
        f >> groupBdrs[g];
        // scale to keV
        groupBdrs[g] *= mec2;
    }

    // Remaining lines: Temperatures and MG data
    // Format:
    // T
    // gto gfrom moment0 [moment1 moment2 ...]
    // <blankline>
    for (UINT iT=0; iT < numTs; ++iT)
    {
        f >> Ts[iT];
        // scale to keV
        Ts[iT] *= mec2;
        const FP T4 = Ts[iT] * Ts[iT] * Ts[iT] * Ts[iT];

        const FP linscale = isnonlin ? nlbase * T4: 1.0;
        const FP renorm = basescale * linscale;

        UINT gfrom = 0U;
        UINT gto = 0U;
        bool finished = false;
        while (!finished)
        {
            // Read one line
            f >> gfrom >> gto;
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
                const UINT loc = gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
                Check(loc < data.size());
                data[loc] = val;
            }

            // Get EOL char
            char c;
            f.get(c);
            Check(c == '\n');

            // Look-ahead to next line to see if a blank line,
            // signalling end of temperature data
            f.get(c);
            if (c == '\n' || c == EOF)
            {
                finished = true;
            }
            else
            {
                f.putback(c);
            }
        }
    }
    f.close();
}

// print contents of struct
void Dense_Compton_Data::print_contents(int verbosity, int precision)
{
    if (verbosity != 0)
        std::cout << "Printing contents at precision " << precision << "...\n";

    // Print sizes
    if (verbosity > 0)
    {
        std::cout << '\n';
        std::cout << "numEvals " << numEvals << '\n';
        std::cout << "numTs " << numTs << '\n';
        std::cout << "numGroups " << numGroups << '\n';
        std::cout << "numLegMoments " << numLegMoments << '\n';
    }

    // Print grids
    if (verbosity > 1)
    {
        std::cout << '\n';
        std::cout << "Group boundaries (keV):\n";
        for (UINT g=0; g < numGroups+1U; ++g)
        {
            if (g > 0)
                std::cout << ' ';
            std::cout << std::setprecision(precision) << groupBdrs[g];
        }
        std::cout << '\n';

        std::cout << "Temperatures (keV):\n";
        for (UINT iT=0; iT < numTs; ++iT)
        {
            if (iT > 0)
                std::cout << ' ';
            std::cout << std::setprecision(precision) << Ts[iT];
        }
        std::cout << '\n';
    }

    // Print contents
    if (verbosity > 2)
    {
        std::cout << '\n';
        std::vector<std::string> evalNames = {"in_lin", "out_lin", "in_nonlin", "out_nonlin", "nldiff"};
        for (UINT eval=0; eval < numEvals; ++eval)
        {
            std::cout << "Eval: " << evalNames[eval] << '\n';
            for (UINT iL=0; iL < numLegMoments; ++iL)
            {
                std::cout << "Legendre moment: " << iL << '\n';

                for (UINT iT=0; iT < numTs; ++iT)
                {
                    std::cout << "Temperature (keV): " << Ts[iT] << '\n';

                    std::cout << "Data (matrix; cm^2/mole):\n";
                    for (UINT gto=0; gto < numGroups; ++gto)
                    {
                        for (UINT gfrom=0; gfrom < numGroups; ++gfrom)
                        {
                            if (gfrom > 0)
                                std::cout << ' ';
                            const UINT loc = gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
                            std::cout << std::setprecision(precision) << data[loc];
                        }
                        std::cout << '\n';
                    }

                    std::cout << "Derivative in T (matrix; cm^2/mole-keV):\n";
                    for (UINT gto=0; gto < numGroups; ++gto)
                    {
                        for (UINT gfrom=0; gfrom < numGroups; ++gfrom)
                        {
                            if (gfrom > 0)
                                std::cout << ' ';
                            const UINT loc = gto + numGroups * (gfrom + numGroups * (iT + numTs * (iL + numLegMoments * eval)));
                            std::cout << std::setprecision(precision) << derivs[loc];
                        }
                        std::cout << '\n';
                    }
                }
            }
        }
    }


    if (verbosity != 0)
    {
        std::cout << '\n';
        std::cout << "...done printing contents\n";
        std::cout << '\n';
    }
}

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
    //std::array<std::string, 1> lins = {"lin"};
    std::array<std::string, 1> lins = {"nonlin"};
#endif

    // Store Compton data in structure
    Dense_Compton_Data dat;

    // Resize
    {
        UINT numfiles = lins.size() * inouts.size();
        std::string const filename = basename + '_' + inouts[0] + '_' + lins[0];
        dat.resize(numfiles, filename);
    }

    // Fill
    UINT counter = 0;
    for (std::string lin : lins) {
        for (std::string inout : inouts) {

            ++counter;
            UINT eval = counter - 1U;

            std::string const filename = basename + '_' + inout + '_' + lin;
            cout << "Reading file: " << filename << endl;

            bool isnonlin = lin.compare("nonlin") == 0;

            dat.read_from_file(eval, filename, isnonlin);
        }
    }

    // Print
    int precision = 3;
    //precision = 16;
    dat.print_contents(verbosity, precision);

    // Check detailed balance
    if (dat.numEvals == 5)
    {
        std::cout << "Detailed balance check...\n";
        if (verbosity <= 0)
            std::cout << "T lindiff/nonlindiff\n";

        // Planck MG integral
        std::vector<FP> bg(dat.numGroups, 0.0);

        // indexes for the evaluations
        // I for inscattering, O for outscattering, f for difference
        // L for linear, N for nonlinear
        std::vector<std::string> evalNames = {"in_lin", "out_lin", "in_nonlin", "out_nonlin", "nldiff"};
        const UINT i_IL = 0;
        const UINT i_OL = 1;
        const UINT i_IN = 2;
        const UINT i_ON = 3;
        //const UINT i_fN = 4;
        // 0th Legendre moment
        const UINT iL = 0;

        if (verbosity > 0)
            std::cout << '\n';
        for (UINT iT = 0; iT < dat.numTs; ++iT)
        {
            // Get T
            const FP T = dat.Ts[iT];
            if (verbosity > 0)
                std::cout << "Temperature (keV): " << std::setprecision(precision) << T << '\n';

            // Compute bg[T]
            if (verbosity > 1)
                std::cout << "Planck spectrum: ";
            FP bgsum = 0.0;
            for (UINT g = 0; g < dat.numGroups; ++g)
            {
                const FP Elow = dat.groupBdrs[g];
                const FP Ehigh = dat.groupBdrs[g+1];
                bg[g] = rtt_cdi::CDI::integratePlanckSpectrum(Elow, Ehigh, T);
                bgsum += bg[g];

                if (verbosity > 1)
                    std::cout << ' ' << std::setprecision(precision) << bg[g];
            }
            if (verbosity > 1)
                std::cout << '\n';
            if (verbosity > 0)
                std::cout << "bgsum: " << std::setprecision(16) << bgsum << '\n';

            // Compute sums for each eval in equilibrium (I=B)
            std::array<FP,4> sums = {0.0, 0.0, 0.0, 0.0};
            for (UINT eval = 0; eval < 4; ++eval)
            {
                for (UINT gfrom = 0; gfrom < dat.numGroups; ++gfrom)
                {
                    FP subsum = 0.0;
                    for (UINT gto = 0; gto < dat.numGroups; ++gto)
                    {
                        // for linear terms, no induced planck[energy_to]
                        const FP bgto = (eval >= 2) ? bg[gto] : 1.0;
                        const FP bgfrom = bg[gfrom];
                        const UINT loc = gto + dat.numGroups * (gfrom + dat.numGroups * (iT + dat.numTs * (iL + dat.numLegMoments * eval)));
                        const FP val = dat.data[loc];
                        subsum += bgto * val * bgfrom;
                    }
                    sums[eval] += subsum;
                }
                if (verbosity > 0)
                    std::cout << evalNames[eval] << " sum: " << sums[eval] << '\n';
            }

            // Print detailed-balance differences
            const FP lindiff = (sums[i_IL] - sums[i_OL]) * 0.75e4 / T;
            const FP nonlindiff = (sums[i_ON] - sums[i_IN]) * 0.75e4 / T;
            const FP ratio = lindiff / nonlindiff;
            if (verbosity > 0)
            {
                std::cout << "lindiff nonlindiff: " << std::setprecision(6) << lindiff << ' ' << nonlindiff << '\n';
                std::cout << "lindiff / nonlindiff: " << ratio << '\n';
            }
            else
                std::cout << T << ' ' << ratio << '\n';

            if (verbosity > 0)
                std::cout << '\n';

        }
        std::cout << "...detailed balance check done\n";
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
        //verbosity = 3;
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
