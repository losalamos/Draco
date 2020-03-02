//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi/CDI_integrateRosselandSpectrum.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "CDI.hh"

namespace rtt_cdi {
//----------------------------------------------------------------------------//
// Rosseland Spectrum Integrators
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
/*!
 * \brief Integrate the Rosseland spectrum over a frequency range.
 *
 * The arguments to this function must all be in consistent units. For example,
 * if low and high are expressed in keV, then the temperature must also be
 * expressed in keV. If low and high are in Hz and temperature is in K, then low
 * and high must first be multiplied by Planck's constant and temperature by
 * Boltzmann's constant before they are passed to this function.
 *
 * \param low lower frequency bound.
 * \param high higher frequency bound.
 * \param T the temperature (must be greater than 0.0)
 *
 * \return integrated normalized Rosseland from low to high
 */
double CDI::integrateRosselandSpectrum(const double low, const double high,
                                       const double T) {
  Require(low >= 0.0);
  Require(high >= low);
  Require(T >= 0.0);

  double planck, rosseland;

  integrate_Rosseland_Planckian_Spectrum(low, high, T, planck, rosseland);

  return rosseland;
}

} // end namespace rtt_cdi

//----------------------------------------------------------------------------//
// end of CDI_integrateRosselandSpectrum.cc
//----------------------------------------------------------------------------//
