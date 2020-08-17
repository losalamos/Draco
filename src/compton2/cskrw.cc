//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   compton2/cskrw.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Basic reader/writer of csk Compton files
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "cdi/CDI.hh" //planck integral
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/XGetopt.hh"
#include "units/PhysicalConstants.hh"
#include "units/UnitSystem.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
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

using UINT = uint64_t;
using FP = double;

struct Dense_Compton_Data {
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
  void compute_nonlinear_difference();
  void compute_temperature_derivatives();
  void write_sparse_binary(std::string fileout);
  void print_contents(int verbosity, int precision);

private:
  void copy_to_sparse(/* ADD */);
  void write_binary(std::string fileout /* ADD */);
};

// resize data and set sizes
// numfiles is number of Compton files; filename is one such file
void Dense_Compton_Data::resize(UINT numfiles, std::string filename) {
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
  groupBdrs.resize(numGroups + 1, 0.0);
  Ts.resize(numTs, 0.0);
  UINT sz = numEvals * numLegMoments * numTs * numGroups * numGroups;
  data.resize(sz, 0.0);
  derivs.resize(sz, 0.0);

  f.close();
}

// read the entire contents of one file
void Dense_Compton_Data::read_from_file(UINT eval, std::string filename,
                                        bool isnonlin) {
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
  const FP nlbase = (4.0 / 9.0) * 2.0 /
                    (hplanck * hplanck * hplanck * cspeed * cspeed) *
                    (mec2 * mec2 * mec2);

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
  for (UINT i = 0; i < numTbreakpoints; ++i) {
    f >> throwaway;
  }

  // Line 3: Group bounds
  for (UINT g = 0; g < numGroups + 1U; ++g) {
    f >> groupBdrs[g];
    // scale to keV
    groupBdrs[g] *= mec2;
  }

  // Remaining lines: Temperatures and MG data
  // Format:
  // T
  // gto gfrom moment0 [moment1 moment2 ...]
  // <blankline>
  for (UINT iT = 0; iT < numTs; ++iT) {
    f >> Ts[iT];
    // scale to keV
    Ts[iT] *= mec2;
    const FP T4 = Ts[iT] * Ts[iT] * Ts[iT] * Ts[iT];

    const FP linscale = isnonlin ? nlbase * T4 : 1.0;
    const FP renorm = basescale * linscale;

    UINT gfrom = 0U;
    UINT gto = 0U;
    bool finished = false;
    while (!finished) {
      // Read one line
      f >> gfrom >> gto;
      // 1-based to 0-based
      gfrom -= 1U;
      gto -= 1U;

      // Read xs
      for (UINT iL = 0; iL < numLegMoments; ++iL) {
        // Read value
        FP val;
        f >> val;
        val *= renorm;

        // Put in 1D data vector
        const UINT loc =
            gto + numGroups *
                      (gfrom +
                       numGroups * (iT + numTs * (iL + numLegMoments * eval)));
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
      if (c == '\n' || c == EOF) {
        finished = true;
      } else {
        f.putback(c);
      }
    }
  }
  f.close();
}

// use 4 evaluations to determine nonlinear difference
// implicit at low E/T; explicit at high E/T
void Dense_Compton_Data::compute_nonlinear_difference() {
  Require(numEvals == 5);

  // Physically, ON - IN.T == IL.T B^-1 - B^-1 OL.
  // Due to quadrature choice and finite-precision arithmetic,
  // the above equality does NOT discretely hold.
  // The LHS suffers from small differences of large numbers at low E/T.
  // The RHS suffers from the same problem when E/T is large and so B is small.
  // To mitigate numerical error and avoid dividing by a B of 0.,
  // we use the RHS at low E/T and use the LHS at high E/T.

  // Determine the cutoff between low and high energies (Ecutoff = N * T)
  // Based on the wgt fxn CSK uses in its MG avg, N <= 25.0
  // Using 25.0 is better for equilibrium stimulated scat
  // Using N>=9.0 is better for non-equil stimulated scat
  const FP Ncutoff = 9.0;

  // Planck MG integral
  std::vector<FP> bg(numGroups, 0.0);

  // indexes for the evaluations
  // I for inscattering, O for outscattering, f for difference
  // L for linear, N for nonlinear
  const UINT e_IL = 0;
  const UINT e_OL = 1;
  const UINT e_IN = 2;
  const UINT e_ON = 3;
  const UINT e_fN = 4;
  // 0th Legendre moment
  const UINT iL = 0;

  for (UINT iT = 0; iT < numTs; ++iT) {
    const FP T = Ts[iT];
    const FP Ecutoff = Ncutoff * T;

    // Compute bg[T]
    FP bgsum = 0.0;
    for (UINT g = 0; g < numGroups; ++g) {
      const FP Elow = groupBdrs[g];
      const FP Ehigh = groupBdrs[g + 1];
      bg[g] = rtt_cdi::CDI::integratePlanckSpectrum(Elow, Ehigh, T);
      bgsum += bg[g];
    }
    // First pass on nldiff
    FP sumlin = 0.0;
    FP sumnonlin = 0.0;
    for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
      for (UINT gto = 0; gto < numGroups; ++gto) {
        // Look at left side of group bounds (use less implicit)
        const FP Eto = groupBdrs[gto];
        const FP Efrom = groupBdrs[gfrom];
        const bool lowE = Eto <= Ecutoff && Efrom <= Ecutoff;

        // Planck is equilibrium distribution
        const FP bgto = bg[gto];
        const FP bgfrom = bg[gfrom];

        std::array<FP, 4> vals;
        // use scattering matrix (no transpose) for outscattering
        for (UINT eval : {e_OL, e_ON}) {
          const UINT loc =
              gto + numGroups *
                        (gfrom + numGroups * (iT + numTs * (iL + numLegMoments *
                                                                     eval)));
          vals[eval] = data[loc];
        }
        // use transpose of scattering matrix for inscattering
        for (UINT eval : {e_IL, e_IN}) {
          const UINT loc =
              gfrom +
              numGroups * (gto + numGroups * (iT + numTs * (iL + numLegMoments *
                                                                     eval)));
          vals[eval] = data[loc];
        }

        // Avoid dividing by zero
        const FP eps = 100 * std::numeric_limits<FP>::min();
        const bool bzero = bgto <= eps || bgfrom <= eps;

        // Take differences of spontaneous and induced rates at equilibrium
        const FP impldiff =
            (bzero) ? 0.0 : vals[e_IL] / bgfrom - vals[e_OL] / bgto;
        const FP expldiff = vals[e_ON] - vals[e_IN];

        // For low E/T, store impldiff; for high E/T, store expldiff
        const UINT loc_fN =
            gto + numGroups *
                      (gfrom +
                       numGroups * (iT + numTs * (iL + numLegMoments * e_fN)));
        data[loc_fN] = (lowE) ? impldiff : expldiff;

        // Keep sums of rates for later ratio
        sumlin += bgto * impldiff * bgfrom;
        sumnonlin += bgto * data[loc_fN] * bgfrom;
      }
    }

    // TODO: Add in `a (c?) T^4` rescaling for fN, ON, and IN

    // Rescale the nonlin diff to get exact detailed balance at equilibrium
    const FP scalenl = sumlin / sumnonlin;
    // we hope scalenl is within a percent or less of 1
    Check(scalenl < 1.2 && scalenl > 0.8);
    for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
      for (UINT gto = 0; gto < numGroups; ++gto) {
        const UINT loc_fN =
            gto + numGroups *
                      (gfrom +
                       numGroups * (iT + numTs * (iL + numLegMoments * e_fN)));
        data[loc_fN] *= scalenl;
      }
    }
  }
}

// Compute and limit temperature derivatives of data
void Dense_Compton_Data::compute_temperature_derivatives() {
  // Check that temperature grid is valid (part 1/2)
  if (numTs < 2) {
    std::fill(derivs.begin(), derivs.end(), 0.0);
    std::cerr << "WARNING: Cannot construct derivatives with only one "
                 "temperature. Aborting routine.";
    std::cerr << std::endl;
    return;
  }

  // Check that temperature grid is valid (part 2/2)
  bool validTs = true;
  for (UINT it = 0; it < (numTs - 1U); ++it) {
    FP T1 = Ts[it];
    FP T2 = Ts[it + 1];
    validTs = validTs && (T1 < T2);
  }
  Insist(validTs,
         "Temperatures are not monotonically increasing and unique.\n");

  // Temporary array for finite differences
  const UINT fd_sz = numGroups * numGroups * (numTs - 1U);
  std::vector<FP> finite_diffs(fd_sz, 0.0);

  for (UINT eval = 0; eval < numEvals; ++eval) {
    for (UINT iL = 0; iL < numLegMoments; ++iL) {
      std::fill(finite_diffs.begin(), finite_diffs.end(), 0.0);

      // Step 1: Compute finite difference in each temperature interval
      for (UINT iT = 0; iT < (numTs - 1U); ++iT) {
        const FP inv_dT = 1.0 / (Ts[iT + 1U] - Ts[iT]);
        for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT gto = 0; gto < numGroups; ++gto) {

            const UINT loc_m =
                gto +
                numGroups *
                    (gfrom +
                     numGroups * (iT + numTs * (iL + numLegMoments * eval)));
            const UINT loc_p =
                gto +
                numGroups *
                    (gfrom + numGroups * ((iT + 1U) +
                                          numTs * (iL + numLegMoments * eval)));
            const FP fd = (data[loc_p] - data[loc_m]) * inv_dT;

            const UINT loc_fd = gto + numGroups * (gfrom + numGroups * iT);
            finite_diffs[loc_fd] = fd;
          }
        }
      }

      // Step 2: Compute derivatives at each temperature point with a
      // limiter based on the finite differences in surrounding intervals

      // Step 2a: First temperature (no limiter; first-order estimate)
      {
        UINT iT = 0;
        const UINT offset =
            numGroups * (iT + numTs * (iL + numLegMoments * eval));
        const UINT offset_fd = numGroups * iT;
        for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT gto = 0; gto < numGroups; ++gto) {
            const UINT loc = gto + numGroups * (gfrom + offset);
            const UINT loc_fd = gto + numGroups * (gfrom + offset_fd);
            derivs[loc] = finite_diffs[loc_fd];
          }
        }
      }

      // Step 2b: Last temperature (no limiter; first-order estimate)
      {
        UINT iT = numTs - 1U;
        const UINT offset =
            numGroups * (iT + numTs * (iL + numLegMoments * eval));
        const UINT offset_fd = numGroups * (iT - 1U);
        for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT gto = 0; gto < numGroups; ++gto) {
            const UINT loc = gto + numGroups * (gfrom + offset);
            const UINT loc_fd = gto + numGroups * (gfrom + offset_fd);
            derivs[loc] = finite_diffs[loc_fd];
          }
        }
      }

      // Step 2c: Interior temperatures (with limiters)
      // m is left of iT; p is right of iT; p and iT have same indices
      for (UINT iT = 1U; iT < (numTs - 1U); ++iT) {
        const FP dT_m = Ts[iT] - Ts[iT - 1U];
        const FP dT_p = Ts[iT + 1U] - Ts[iT];
        const FP f_m = (2. * dT_m + dT_p) / (3. * (dT_m + dT_p));
        const FP f_p = 1.0 - f_m;

        const UINT offset =
            numGroups * (iT + numTs * (iL + numLegMoments * eval));

        const UINT offset_fd_m = numGroups * (iT - 1U);
        const UINT offset_fd_p = numGroups * iT;

        for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
          for (UINT gto = 0; gto < numGroups; ++gto) {
            const UINT loc_fd_m = gto + numGroups * (gfrom + offset_fd_m);
            const UINT loc_fd_p = gto + numGroups * (gfrom + offset_fd_p);
            const FP fd_m = finite_diffs[loc_fd_m];
            const FP fd_p = finite_diffs[loc_fd_p];

            // For more information on method, see
            // scipy.interpolate.pchip and
            // http://dx.doi.org/10.1137/1.9780898717952

            const int sign_m = int(0.0 < fd_m) - int(fd_m < 0.0);
            const int sign_p = int(0.0 < fd_p) - int(fd_p < 0.0);

            const FP d = (sign_m * sign_p > 0)
                             ? (fd_m * fd_p) / (f_m * fd_m + f_p * fd_p)
                             : 0.0;
            const UINT loc = gto + numGroups * (gfrom + offset);
            derivs[loc] = d;
          }
        }
      }
    }
  }
}

// Sparsify data and print to binary
void Dense_Compton_Data::write_sparse_binary(std::string fileout) {
  std::cout << "Sparse and binary not yet implemented.\n";
  copy_to_sparse();
  write_binary(fileout);
}

// Sparsify data
void Dense_Compton_Data::copy_to_sparse() {
  // RM
  std::cout << "  to sparse\n";
}

// Write to binary
void Dense_Compton_Data::write_binary(std::string fileout) {
  // RM
  std::cout << "  write binary\n";

  auto fout = std::ofstream(fileout, std::ios::out | std::ios::binary);

  // binary type
  char filetype[] = " csk ";
  fout.write(&filetype[0], sizeof(filetype));

  fout.close();
}

// print contents of struct
void Dense_Compton_Data::print_contents(int verbosity, int precision) {
  if (verbosity != 0)
    std::cout << "Printing contents at precision " << precision << "...\n";

  // Print sizes
  if (verbosity > 0) {
    std::cout << '\n';
    std::cout << "numEvals " << numEvals << '\n';
    std::cout << "numTs " << numTs << '\n';
    std::cout << "numGroups " << numGroups << '\n';
    std::cout << "numLegMoments " << numLegMoments << '\n';
  }

  // Print grids
  if (verbosity > 1) {
    std::cout << '\n';
    std::cout << "Group boundaries (keV):\n";
    for (UINT g = 0; g < numGroups + 1U; ++g) {
      if (g > 0)
        std::cout << ' ';
      std::cout << std::setprecision(precision) << groupBdrs[g];
    }
    std::cout << '\n';

    std::cout << "Temperatures (keV):\n";
    for (UINT iT = 0; iT < numTs; ++iT) {
      if (iT > 0)
        std::cout << ' ';
      std::cout << std::setprecision(precision) << Ts[iT];
    }
    std::cout << '\n';
  }

  // Print contents
  if (verbosity > 2) {
    std::cout << '\n';
    std::vector<std::string> evalNames = {"in_lin", "out_lin", "in_nonlin",
                                          "out_nonlin", "nldiff"};
    for (UINT eval = 0; eval < numEvals; ++eval) {
      std::cout << "Eval: " << evalNames[eval] << '\n';
      for (UINT iL = 0; iL < numLegMoments; ++iL) {
        std::cout << "Legendre moment: " << iL << '\n';

        for (UINT iT = 0; iT < numTs; ++iT) {
          std::cout << "Temperature (keV): " << Ts[iT] << '\n';

          std::cout << "Data (matrix; cm^2/mole):\n";
          for (UINT gto = 0; gto < numGroups; ++gto) {
            if (gto > 0)
              continue; // HACK!
            for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
              if (gfrom > 0)
                std::cout << ' ';
              const UINT loc =
                  gto +
                  numGroups *
                      (gfrom +
                       numGroups * (iT + numTs * (iL + numLegMoments * eval)));
              std::cout << std::setprecision(precision) << data[loc];
            }
            std::cout << '\n';
          }

          std::cout << "Derivative in T (matrix; cm^2/mole-keV):\n";
          for (UINT gto = 0; gto < numGroups; ++gto) {
            if (gto > 0)
              continue; // HACK!
            for (UINT gfrom = 0; gfrom < numGroups; ++gfrom) {
              if (gfrom > 0)
                std::cout << ' ';
              const UINT loc =
                  gto +
                  numGroups *
                      (gfrom +
                       numGroups * (iT + numTs * (iL + numLegMoments * eval)));
              std::cout << std::setprecision(precision) << derivs[loc];
            }
            std::cout << '\n';
          }
        }
      }
    }
  }

  if (verbosity != 0) {
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

#if 0
    std::array<std::string, 2> inouts = {"in", "out"};
    std::array<std::string, 2> lins = {"lin", "nonlin"};
#else
  std::array<std::string, 1> inouts = {"in"};
  std::array<std::string, 1> lins = {"lin"};
  //std::array<std::string, 1> lins = {"nonlin"};
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

  // Use filled data to compute derived data
  // (4 evals from raw data plus 1 eval of derived data)
  if (dat.numEvals == 5) {
    dat.compute_nonlinear_difference();
  }
  dat.compute_temperature_derivatives();

  // Save to binary
  std::string fileout = basename + "_b";
  std::cout << "Writing file: " << fileout << '\n';
  dat.write_sparse_binary(fileout);

  // Print
  int precision = 3;
  //precision = 16;
  dat.print_contents(verbosity, precision);

  // RM!
  return;

  // Check detailed balance
  if (dat.numEvals == 5) {
    std::cout << "Detailed balance check...\n";
    if (verbosity <= 0)
      std::cout << "T lindiff/nonlindiff\n";

    // Planck MG integral
    std::vector<FP> bg(dat.numGroups, 0.0);

    // indexes for the evaluations
    // I for inscattering, O for outscattering, f for difference
    // L for linear, N for nonlinear
    std::vector<std::string> evalNames = {"in_lin", "out_lin", "in_nonlin",
                                          "out_nonlin", "nldiff"};
    const UINT e_IL = 0;
    const UINT e_OL = 1;
    const UINT e_IN = 2;
    const UINT e_ON = 3;
    const UINT e_fN = 4;
    // 0th Legendre moment
    const UINT iL = 0;

    if (verbosity > 0)
      std::cout << '\n';
    for (UINT iT = 0; iT < dat.numTs; ++iT) {
      // Get T
      const FP T = dat.Ts[iT];
      if (verbosity > 0)
        std::cout << "Temperature (keV): " << std::setprecision(precision) << T
                  << '\n';

      // Compute bg[T]
      if (verbosity > 1)
        std::cout << "Planck spectrum: ";
      FP bgsum = 0.0;
      for (UINT g = 0; g < dat.numGroups; ++g) {
        const FP Elow = dat.groupBdrs[g];
        const FP Ehigh = dat.groupBdrs[g + 1];
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
      std::array<FP, 5> sums = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (UINT eval = 0; eval < sums.size(); ++eval) {
        for (UINT gfrom = 0; gfrom < dat.numGroups; ++gfrom) {
          FP subsum = 0.0;
          for (UINT gto = 0; gto < dat.numGroups; ++gto) {
            // for linear terms, no induced planck[energy_to]
            const FP bgto = (eval >= 2) ? bg[gto] : 1.0;
            const FP bgfrom = bg[gfrom];
            const UINT loc =
                gto +
                dat.numGroups *
                    (gfrom +
                     dat.numGroups *
                         (iT + dat.numTs * (iL + dat.numLegMoments * eval)));
            const FP val = dat.data[loc];
            subsum += bgto * val * bgfrom;
          }
          sums[eval] += subsum;
        }
        if (verbosity > 1)
          std::cout << evalNames[eval] << " sum: " << sums[eval] << '\n';
      }

      // Print detailed-balance differences
      const FP scale = 0.75e4 / T;
      const FP lindiff = (sums[e_IL] - sums[e_OL]) * scale;
      const FP nonlindiff_raw = (sums[e_ON] - sums[e_IN]) * scale;
      const FP nonlindiff_use = sums[e_fN] * scale;
      const FP ratio_raw = lindiff / nonlindiff_raw;
      const FP ratio_use = lindiff / nonlindiff_use;
      if (verbosity > 1)
        std::cout << "lindiff nonlindiff-RAW nonlindiff-USE: "
                  << std::setprecision(6) << lindiff << ' ' << nonlindiff_raw
                  << ' ' << nonlindiff_use << '\n';
      if (verbosity > 0) {
        std::cout << "lindiff / nonlindiff_raw: " << ratio_raw << '\n';
        std::cout << "lindiff / nonlindiff_use: " << ratio_use << '\n';
      } else
        std::cout << T << ' ' << ratio_use << '\n';

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
    verbosity = 2;
    verbosity = 3;
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
