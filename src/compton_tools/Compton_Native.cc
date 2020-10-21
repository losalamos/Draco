//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   compton_tools/Compton_Native.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Implementation file for native compton bindary-read and temperature interpolation
 * \note   Copyright (C) 2017-2020 Triad National Security, LLC. All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "compton_tools/Compton_Native.hh"
#include "c4/C4_Functions.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <iostream>

using UINT = size_t;
using FP = double;
using vec = std::vector<FP>;

namespace rtt_compton_tools {

//------------------------------------------------------------------------------------------------//
// FREE-FLOATING HELPER FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief find location in a sorted list and min/max value to be within list
 *
 * \param[in] xs vector of monotonically increasing and unique values
 * \param[inout] x value whose location in xs needs to be found;
 *               modified so that xs[index] <= x <= xs[index+1]
 * \return Index in [0, xs.size()-2] so that xs[index] and xs[index+1] are valid
 *
 * Use a binary search to find the location of x in xs
 */
UINT find_index(const std::vector<FP> &xs, FP &x) {
  Require(xs.size() >= 2);

  // Do binary search for index
  auto loc = std::upper_bound(xs.begin(), xs.end(), x);
  // loc - begin guaranteed to be in [0,len]
  UINT index = (loc - xs.cbegin());

  // Move index and x to be interior
  UINT len = xs.size();
  FP xmin = xs[0];
  FP xmax = xs[len - 1U];
  index = std::min(len - 1U, std::max(static_cast<UINT>(1), index)) - 1U;
  x = std::max(xmin, std::min(xmax, x));
  return index;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief compute hermite polynomial for a given value and the left/right grid points
 *
 * \param[in] x value of independent variable
 * \param[in] xL gridpoint to the left of x (xL <= x)
 * \param[in] xR gridpoint to the right of x (x <= xR)
 * \return length-4 array of Hermite polynomials used for interpolation from data at xL and xR to x
 *
 * If vL/vR are values of a function at xL/xR and dL/dR are derivatives of a function at xL/xR, then
 * H[0] * vL + H[1] * vR + H[2] * dL + H[3] * dR
 * interpolates the value of the function at x, where this function computes and returns H
 */
std::array<FP, 4> hermite(FP x, FP xL, FP xR) {
  Require(xL <= x && x <= xR && xL < xR);

  // Spacing of interval
  FP dx = xR - xL;
  // Left/right basis functions
  FP bL = (xR - x) / dx;
  FP bR = (x - xL) / dx;

  // Hermite functions
  std::array<FP, 4> hermite_eval{bL * bL * (3 - 2 * bL), bR * bR * (3 - 2 * bR),
                                 -dx * bL * bL * (bL - 1), dx * bR * bR * (bR - 1)};
  return hermite_eval;
}

//------------------------------------------------------------------------------------------------//
// SETUP FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Constructor for Compton_Native
 *
 * \param[in] filename Name of the binary csk data file to be read
 *
 * This function has rank 0 read the binary csk data file, which fills in the class' data members.
 * The data members are then broadcast to other MPI ranks to finish their construction.
 */
Compton_Native::Compton_Native(const std::string &filename)
    : num_temperatures_(0U), num_groups_(0U), num_leg_moments_(0U), num_evals_(0U), num_points_(0U),
      Ts_(0U), Egs_(0U), first_groups_(0U), indexes_(0U), data_(0U), derivs_(0U) {
  Require(filename.length() > 0U);
  int rank = rtt_c4::node();
  constexpr int bcast_rank = 0;
  int errcode = 0;

  if (rank == bcast_rank) {
    errcode = read_binary(filename);
  }
  broadcast_MPI(errcode);

  Ensure(check_class_invariants());
}

//------------------------------------------------------------------------------------------------//
// BINARY READING AND BROADCASTING FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Helper member function that broadcasts csk data from rank 0 to all ranks and sets data
 *
 * \param[in] errcode If non-zero, everyone aborts.
 *
 * Uses rtt_c4's broadcast to send arrays and vectors to all ranks with the sparse csk data,
 * temperature / energy grids, and sizes
 */
void Compton_Native::broadcast_MPI(int errcode) {
  rtt_c4::global_max(errcode);
  Insist(errcode == 0, "Non-zero errorcode. Exiting.");

  int rank = rtt_c4::node();
  constexpr int bcast_rank = 0;

  // Broadcast sizes
  UINT data_size = data_.size();
  std::array<UINT, 6> pack = {num_temperatures_, num_groups_, num_leg_moments_,
                              num_evals_,        num_points_, data_size};
  rtt_c4::broadcast(&pack[0], pack.size(), bcast_rank);
  UINT p = 0;
  num_temperatures_ = pack[p++];
  num_groups_ = pack[p++];
  num_leg_moments_ = pack[p++];
  num_evals_ = pack[p++];
  num_points_ = pack[p++];
  data_size = pack[p++];

  // Derived sizes
  const auto tsz = static_cast<int>(num_temperatures_);
  const auto egsz = static_cast<int>(num_groups_ + 1U);
  const auto fgsz = static_cast<int>(num_temperatures_ * num_groups_);
  const auto isz = static_cast<int>(fgsz + 1U);
  const auto dsz = static_cast<int>(data_size);

  // Broadcast grids
  if (rank != bcast_rank)
    Ts_.resize(tsz);
  rtt_c4::broadcast(&Ts_[0], tsz, 0);
  if (rank != bcast_rank)
    Egs_.resize(egsz);
  rtt_c4::broadcast(&Egs_[0], egsz, 0);

  // Broadcast sparse data structures
  if (rank != bcast_rank)
    first_groups_.resize(fgsz);
  rtt_c4::broadcast(&first_groups_[0], fgsz, 0);
  if (rank != bcast_rank)
    indexes_.resize(isz);
  rtt_c4::broadcast(&indexes_[0], isz, 0);

  // Broadcast data itself
  if (rank != bcast_rank)
    data_.resize(dsz);
  rtt_c4::broadcast(&data_[0], dsz, 0);
  if (rank != bcast_rank)
    derivs_.resize(dsz);
  rtt_c4::broadcast(&derivs_[0], dsz, 0);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Helper member function to read a binary csk file and set class data
 *
 * \param[in] filename Path to csk binary file
 * \return errcode Zero if read is successful, otherwise non-zero
 *
 * Reads a binary csk file by interpreting the characters as 64-bit unsigned ints and doubles
 */
int Compton_Native::read_binary(const std::string &filename) {

  // Read
  auto fin = std::ifstream(filename, std::ios::in | std::ios::binary);

  // Check for valid file stream
  try {
    fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  } catch (std::fstream::failure &e) {
    return 2;
  }

  // Ensure valid type
  char expected[] = " csk ";
  char file_type[sizeof(expected)];
  fin.read(file_type, sizeof(file_type));
  if (std::strcmp(file_type, expected) != 0) {
    std::cerr << "Expecting binary file " << filename << " to start with '" << expected
              << "' but got '" << file_type << "'";
    std::cerr << std::endl;
    return 1;
  }

  UINT binary_ordering;
  UINT version_major;
  UINT version_minor;
  fin.read(reinterpret_cast<char *>(&version_major), sizeof(UINT));
  fin.read(reinterpret_cast<char *>(&version_minor), sizeof(UINT));
  fin.read(reinterpret_cast<char *>(&binary_ordering), sizeof(UINT));
  if (version_major != 1 || binary_ordering > 1) {
    std::cerr << "Expecting a CSK binary file (version 1) with ordering 0 or 1 "
                 "but got "
              << version_major << " with ordering " << binary_ordering;
    std::cerr << std::endl;
    return 2;
  }

  constexpr UINT n = 7;
  std::array<UINT, n> szs;
  for (UINT i = 0; i < n; ++i)
    fin.read(reinterpret_cast<char *>(&szs[i]), sizeof(szs[i]));
  UINT j = 0;
  UINT tsz = szs[j++];
  UINT gsz = szs[j++];
  UINT lsz = szs[j++];
  UINT esz = szs[j++];
  UINT fgsz = szs[j++];
  UINT isz = szs[j++];
  UINT dsz = szs[j++];

  num_temperatures_ = tsz;
  num_groups_ = gsz;
  num_leg_moments_ = lsz;
  num_evals_ = esz;
  // point is (leg moment, eval) pair
  // first eval has all leg moments; others have only the 0th moment
  num_points_ = num_evals_ + num_leg_moments_ - 1U;
  UINT egsz = gsz + 1;

#if 0
    std::cout << "DBG num_temperatures_ " << num_temperatures_ << '\n';
    std::cout << "DBG num_groups_ " << num_groups_ << '\n';
    std::cout << "DBG num_leg_moments_ " << num_leg_moments_ << '\n';
    std::cout << "DBG num_evals_ " << num_evals_ << '\n';
    std::cout << "DBG num_points_ " << num_points_ << '\n';
    std::cout << "DBG len(first_groups_) " << fgsz << '\n';
    std::cout << "DBG len(indexes_) " << isz << '\n';
    std::cout << "DBG len(data_/derivs_) " << dsz << '\n';
#endif

  Ts_.resize(tsz);
  for (UINT i = 0; i < tsz; ++i)
    fin.read(reinterpret_cast<char *>(&Ts_[i]), sizeof(Ts_[i]));

  Egs_.resize(egsz);
  for (UINT i = 0; i < egsz; ++i)
    fin.read(reinterpret_cast<char *>(&Egs_[i]), sizeof(Egs_[i]));

  first_groups_.resize(fgsz);
  for (UINT i = 0; i < fgsz; ++i)
    fin.read(reinterpret_cast<char *>(&first_groups_[i]), sizeof(first_groups_[i]));

  indexes_.resize(isz);
  for (UINT i = 0; i < isz; ++i)
    fin.read(reinterpret_cast<char *>(&indexes_[i]), sizeof(indexes_[i]));

  data_.resize(dsz);
  for (UINT i = 0; i < dsz; ++i)
    fin.read(reinterpret_cast<char *>(&data_[i]), sizeof(data_[i]));

  derivs_.resize(dsz);
  for (UINT i = 0; i < dsz; ++i)
    fin.read(reinterpret_cast<char *>(&derivs_[i]), sizeof(derivs_[i]));

  fin.close();

#if 0
    std::cout << '\n';
    std::cout << "-------------------------------------------------\n";

    std::cout << "\nDBG Ts";
    for (UINT i = 0; i < tsz; ++i)
      std::cout << ", " << Ts_[i] / 510.998;
    std::cout << '\n';

    std::cout << "\nDBG Egs";
    for (UINT i = 0; i < egsz; ++i)
      std::cout << ", " << Egs_[i] / 510.998;
    std::cout << '\n';

    std::cout << "\nDBG first_groups";
    for (UINT i = 0; i < fgsz; ++i)
      std::cout << ", " << first_groups_[i];
    std::cout << '\n';

    std::cout << "\nDBG indexes";
    for (UINT i = 0; i < isz; ++i)
      std::cout << ", " << indexes_[i];
    std::cout << '\n';

    std::cout << "\nDBG data";
    for (UINT p = 0; p < num_points_; ++p) {
      std::cout << "\nPoint " << p << "\n";
      UINT szp = dsz / num_points_;
      for (UINT ii = 0; ii < fgsz; ++ii) {
        UINT istrt = indexes_[ii] + p * szp;
        UINT iend = indexes_[ii + 1] + p * szp;
        std::cout << "  index " << ii;
        for (UINT i = istrt; i < iend; ++i) {
          std::cout << std::setprecision(12); // << std::scientific;
          std::cout << ", " << data_[i] / 0.075116337052433;
        }
        std::cout << '\n';
      }
    }
    std::cout << std::defaultfloat;

    std::cout << "-------------------------------------------------\n";
    std::cout << '\n';
#endif

  return 0;
}

//------------------------------------------------------------------------------------------------//
// TEMPERATURE INTERPOLATION FUNCTIONS
//------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate csk data in temperature and return dense linear inscattering matrix
 *
 * \param[inout] inscat The flattened (1D), dense inscattering matrix as a vector
 *               Order of inscat is (slow) [moment, group-to, group-from] (fast)
 *               Does NOT need to be the right size prior to calling
 * \param[in] Te_keV The electron temperature in keV at which the interpolation is desired
 * \param[in] num_moments_truncate The maximum number of Legendre moments to use
 *            Function will use the minimum of this variable and number of moments in the data
 */
void Compton_Native::interp_dense_inscat(vec &inscat, double Te_keV,
                                         size_t num_moments_truncate) const {
  // Finds index and nudges Teff such that
  // Ts_[index] <= Teff <= Ts_[index+1] and 0 <= index <= Ts_.size()-2;
  FP Teff = Te_keV;
  UINT iT = rtt_compton_tools::find_index(Ts_, Teff);

  // Fill Hermite function
  std::array<FP, 4> hermite = rtt_compton_tools::hermite(Teff, Ts_[iT], Ts_[iT + 1U]);

  // Precompute some sparse indexes
  const UINT sz = indexes_[indexes_.size() - 1];
  const UINT end_leg = std::min(static_cast<UINT>(num_moments_truncate), num_leg_moments_);
  const UINT eval_offset = 0; // in_lin

  // Resize and fill with zeros
  inscat.resize(end_leg * num_groups_ * num_groups_);
  std::fill(inscat.begin(), inscat.end(), 0.0);

  // Apply Hermite function
  for (UINT k = 0; k < end_leg; ++k) {
    for (UINT gfrom = 0; gfrom < num_groups_; ++gfrom) {
      const UINT offset_jj = gfrom + num_groups_ * num_groups_ * k;
      // Get contributions from both Ts_[iT] and Ts_[iT+1]
      for (UINT n = 0; n < 2U; ++n) {
        const UINT i = gfrom + num_groups_ * (iT + n);
        const UINT first_gto = first_groups_[i];
        const UINT num_entries = indexes_[i + 1U] - indexes_[i];
        const UINT offset_ii = indexes_[i] + sz * k + eval_offset;
        for (UINT dg = 0; dg < num_entries; ++dg) {
          const UINT gto = dg + first_gto;
          const UINT ii = dg + offset_ii;
          const UINT jj = gto * num_groups_ + offset_jj;
          inscat[jj] += hermite[0U + n] * data_[ii] + hermite[2U + n] * derivs_[ii];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate csk data in temperature and return linear outscattering vector
 *
 * \param[inout] outscat The 1D linear outscattering array at the desired temperature
 *               Has been summed over outgoing group so only index is [group-from]
 *               Does NOT need to be the right size prior to calling
 * \param[in] Te_keV The electron temperature in keV at which the interpolation is desired
 */
void Compton_Native::interp_linear_outscat(vec &outscat, double Te_keV) const {
  // Finds index and nudges Teff st Ts_[index] <= Teff <= Ts_[index+1]
  // and 0 <= index <= Ts_.size()-2;
  FP Teff = Te_keV;
  UINT iT = rtt_compton_tools::find_index(Ts_, Teff);

  // Fill Hermite function
  std::array<FP, 4> hermite = rtt_compton_tools::hermite(Teff, Ts_[iT], Ts_[iT + 1U]);

  // Precompute some sparse indexes
  const UINT sz = indexes_[indexes_.size() - 1];
  const UINT eval_offset = sz * num_leg_moments_; // out_lin

  // Resize and fill with zeros
  outscat.resize(num_groups_);
  std::fill(outscat.begin(), outscat.end(), 0.0);

  // Apply Hermite function
  for (UINT gfrom = 0; gfrom < num_groups_; ++gfrom) {
    // Get contributions from both Ts_[iT] and Ts_[iT+1U]
    for (UINT n = 0; n < 2U; ++n) {
      const UINT i = gfrom + num_groups_ * (iT + n);
      const UINT num_entries = indexes_[i + 1U] - indexes_[i];
      const UINT offset = indexes_[i] + eval_offset;
      for (UINT dg = 0; dg < num_entries; ++dg) {
        const UINT ii = dg + offset;
        outscat[gfrom] += hermite[0U + n] * data_[ii] + hermite[2U + n] * derivs_[ii];
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Interpolate csk data in temperature and add nonlinear difference to outscattering
 *
 * \param[inout] outscat The 1D net outscattering array at the desired temperature
 *               Has been summed over outgoing group so only index is [group-from]
 *               MUST be the right size (# groups) and initialized with data prior to calling
 * \param[in] Te_keV The electron temperature in keV at which the interpolation is desired
 * \param[in] phi The multigroup radiation field of size number of groups
 * \param[in] scale The scale for phi:
 *                  when the radiation is in equilibrium, sum_g phi_g = scale
 *
 * Adds the difference (nonlinear outscattering minus nonlinear inscattering) to the outscattering
 * vector. The contribution is nonlinear because it depends on phi, the radiation field.
 * The use of scale allows phi to be passed in with arbitrary normalization (4pi, c, a, etc.).
 */
void Compton_Native::interp_nonlin_diff_and_add(vec &outscat, double Te_keV,
                                                const std::vector<double> &phi,
                                                double scale) const {
  // Adds to existing outscat vector
  Require(outscat.size() == num_groups_);
  Require(phi.size() == num_groups_);

  // Finds index and nudges Teff st Ts_[index] <= Teff <= Ts_[index+1]
  // and 0 <= index <= Ts_.size()-2;
  FP Teff = Te_keV;
  UINT iT = rtt_compton_tools::find_index(Ts_, Teff);

  // Fill Hermite function
  std::array<FP, 4> hermite = rtt_compton_tools::hermite(Teff, Ts_[iT], Ts_[iT + 1U]);

  // Precompute some sparse indexes
  const UINT sz = indexes_[indexes_.size() - 1];
  const UINT eval_offset = sz * (num_leg_moments_ + 1U); // nl_diff

  // Precompute constants
  const FP invscale = scale > 0.0 ? static_cast<FP>(1.0 / scale) : 0.0;

  // Apply Hermite function
  for (UINT gfrom = 0; gfrom < num_groups_; ++gfrom) {
    // Get contributions from both Ts_[iT] and Ts_[iT+1U]
    for (UINT n = 0; n < 2U; ++n) {
      const UINT i = gfrom + num_groups_ * (iT + n);
      const UINT first_gto = first_groups_[i];
      const UINT num_entries = indexes_[i + 1U] - indexes_[i];
      const UINT offset = indexes_[i] + eval_offset;
      for (UINT dg = 0; dg < num_entries; ++dg) {
        const UINT gto = dg + first_gto;
        const FP mag = invscale * phi[gto];
        const UINT ii = dg + offset;
        const FP val = hermite[0U + n] * data_[ii] + hermite[2U + n] * derivs_[ii];
        outscat[gfrom] += mag * val;
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
// UNIMPLEMENTED, POTENTIAL FUTURE FUNCTIONS
//------------------------------------------------------------------------------------------------//

#if 0
void Compton_Native::interp_matvec(vec &x, const vec &leftscale,
                             const vec &rightscale, double Te_keV,
                             bool zeroth_moment_only) const {
  // TODO: Redo interface? Not in-place??

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_matvec_transpose(vec &xT, const vec &leftscale,
                                       const vec &rightscale, double Te_keV,
                                       bool zeroth_moment_only) const {
  // TODO: Redo interface? Not in-place??

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_sparse_inscat(Sparse_Compton_Matrix &inscat,
                                    const vec &leftscale, const vec &rightscale,
                                    double Te_keV,
                                    bool zeroth_moment_only) const {
  // TODO: Allow offset to fill in data?

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_dense_inscat(vec &inscat, const vec &leftscale,
                                   const vec &rightscale, double Te_keV,
                                   bool zeroth_moment_only) const {
  // Ordering of inscat is 1D array (slow) [moment, group-to, group-from] (fast)
  const vec &L = leftscale;
  const vec &R = rightscale;
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;

  inscat.resize(num_groups_ * num_groups_ * end_leg);
  std::fill(inscat.begin(), inscat.end(), 0.0);

  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_linear_outscat(vec &outscat, const vec &leftscale,
                                     const vec &rightscale,
                                     double Te_keV) const {
  outscat.resize(num_groups_);
  std::fill(outscat.begin(), outscat.end(), 0.0);

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  // TODO: implement
}

//------------------------------------------------------------------------------------------------//

void Compton_Native::interp_nonlinear_diff(vec &nldiff, const vec &leftscale,
                                     const vec &rightscale, const vec &flux,
                                     double flux_scale, double Te_keV) const {
  // Need "const" that changes based on DBC level?
  // Require(leftscale.len() == num_groups_);
  // Require(rightscale.len() == num_groups_);
  // Require(flux.len() == num_groups_);
  // Require(flux_scale > 0.0);
  // Require(Te_keV >= 0.0);

  nldiff.resize(num_groups_);
  std::fill(nldiff.begin(), nldiff.end(), 0.0);

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const fscale = flux_scale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  // TODO: implement

  // Ensure(nldiff.len() == num_groups_);
}
#endif
//------------------------------------------------------------------------------------------------//

} // namespace rtt_compton_tools

//------------------------------------------------------------------------------------------------//
// End compton_tools/Compton_Native.cc
//------------------------------------------------------------------------------------------------//