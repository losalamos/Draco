//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI_integratePlanckSpectrum.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: CDI.cc 7388 2015-01-22 16:02:07Z kellyt $
//---------------------------------------------------------------------------//

#include "CDI.hh"

namespace rtt_cdi {
//---------------------------------------------------------------------------//
/*!
 * \brief Integrate the Planckian spectrum over a frequency range.
 *
 * The arguments to this function must all be in consistent units. For
 * example, if low and high are expressed in keV, then the temperature must
 * also be expressed in keV. If low and high are in Hz and temperature is in
 * K, then low and high must first be multiplied by Planck's constant and
 * temperature by Boltzmann's constant before they are passed to this function.
 *
 * \param low lower frequency bound.
 * \param high higher frequency bound.
 * \param T the temperature (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian from low to high
 *
 */
double CDI::integratePlanckSpectrum(double low, double high, const double T) {
  Require(low >= 0.0);
  Require(high >= low);
  Require(T >= 0.0);

  // return 0 if temperature is a hard zero
  if (T == 0.0)
    return 0.0;

  // Sale the frequencies by temperature
  low /= T;
  high /= T;

  double integral = integrate_planck(high) - integrate_planck(low);

  Ensure(integral >= 0.0);
  Ensure(integral <= 1.0);

  return integral;
}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
// end of CDI_integratePlanckSpectrum.cc
//---------------------------------------------------------------------------//
