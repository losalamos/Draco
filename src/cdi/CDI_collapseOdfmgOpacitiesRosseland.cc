//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI_collapseOdfmgOpacitiesRosseland.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "CDI.hh"
#include "ds++/Safe_Divide.hh"
#include <numeric>

namespace rtt_cdi {

//---------------------------------------------------------------------------//
/*!
 * \brief Collapse a multigroup-multiband opacity set into a single
 *        representative value weighted by the Rosseland function.
 *
 * \param groupBounds The vector of group boundaries. Size n+1
 * \param opacity   A vector of multigroup opacity data.
 * \param rosselandSpectrum A vector of Rosseland integrals for all groups in
 *                  the spectrum (normally generated via
 *                  CDI::integrate_Rosseland_Planckian_Sectrum(...).
 * \param bandWidths The width of odf-bands.
 * \return A single interval Rosseland weighted opacity value.
 *
 * Typically, CDI::integrate_Rosseland_Planckian_Spectrum is called before this
 * function to obtain rosselandSpectrum.
 *
 * There are 2 special cases that we check for:
 * 1. All opacities are zero - just return 0.0;
 * 2. The Rosseland Integral is very small (or zero).  In this case, perform a
 *    modified calculation that does not depend on the Rosseland integral.
 *
 * If neither of the special cases are in effect, then do the normal evaluation.
 */
double CDI::collapseOdfmgOpacitiesRosseland(
    std::vector<double> const &groupBounds,
    std::vector<std::vector<double>> const &opacity,
    std::vector<double> const &rosselandSpectrum,
    std::vector<double> const &bandWidths) {
  Require(groupBounds.size() > 0);
  Require(opacity.size() == groupBounds.size() - 1);
  Require(opacity[0].size() == bandWidths.size());
  Require(rosselandSpectrum.size() == groupBounds.size() - 1);

  // If all opacities are zero, then the Rosseland mean will also be zero.
  double const eps(std::numeric_limits<double>::epsilon());
  double opacity_sum(0.0);
  for (size_t g = 1; g < groupBounds.size(); ++g)
    opacity_sum +=
        std::accumulate(opacity[g - 1].begin(), opacity[g - 1].end(), 0.0);
  if (rtt_dsxx::soft_equiv(opacity_sum, 0.0, eps)) {
    // std::cerr << "\nWARNING (CDI.cc::"
    //           << "collapseMultigroupOpacitiesRosseland)::"
    //           << "\n\tComputing Rosseland Opacity when all opacities"
    //           << " are zero!" << std::endl;
    return 0.0;
  }

  // Integrate the unnormalized Rosseland over the group spectrum
  // int_{\nu_0}^{\nu_G}{d\nu dB(\nu,T)/dT}
  double const rosseland_integral =
      std::accumulate(rosselandSpectrum.begin(), rosselandSpectrum.end(), 0.0);

  // If the group bounds are well outside the Rosseland Spectrum at the current
  // temperature, our algorithm may return a value that is within machine
  // precision of zero.  In this case, we assume that this occurs when the
  // temperature -> 0, so that limit(T->0) dB/dT = \delta(\nu).  In this case we
  // have:
  //
  // sigma_R = sigma(g=0)

  // Initialize sum
  double inv_sig_r_sum(0.0);

  size_t const numGroups = groupBounds.size() - 1;
  size_t const numBands = bandWidths.size();

  if (rosseland_integral < eps) {
    for (size_t ib = 1; ib <= numBands; ++ib) {
      Check(opacity[0][ib - 1] >= 0.0);
      inv_sig_r_sum +=
          rtt_dsxx::safe_pos_divide(bandWidths[ib - 1], opacity[0][ib - 1]);
    }
    return 1.0 / inv_sig_r_sum;
  }
  Check(rosseland_integral > 0.0);

  // Perform integration of (1/sigma) * d(b_g)/dT over all groups:
  // int_{\nu_0}^{\nu_G}{d\nu (1/sigma(\nu,T)) * dB(\nu,T)/dT}

  // Rosseland opacity:

  //    1      int_{\nu_0}^{\nu_G}{d\nu (1/sigma(\nu,T)) * dB(\nu,T)/dT}
  // ------- = ----------------------------------------------------------
  // sigma_R   int_{\nu_0}^{\nu_G}{d\nu dB(\nu,T)/dT}

  // Accumulated quantities for the Rosseland opacities:
  for (size_t g = 1; g <= numGroups; ++g) {
    Check(rosselandSpectrum[g - 1] >= 0.0);
    for (size_t ib = 1; ib <= numBands; ++ib) {
      Check(opacity[g - 1][ib - 1] >= 0.0);
      Check((g - 1) * numBands + ib - 1 < numBands * numGroups);

      inv_sig_r_sum += rtt_dsxx::safe_pos_divide(rosselandSpectrum[g - 1] *
                                                     bandWidths[ib - 1],
                                                 opacity[g - 1][ib - 1]);
    }
  }
  Check(inv_sig_r_sum > 0.0);
  return rosseland_integral / inv_sig_r_sum;
}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
// end of CDI_collapseOdfmgOpacitiesRosseland.cc
//---------------------------------------------------------------------------//
