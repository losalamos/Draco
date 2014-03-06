//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI_integrateRosselandSpectrum.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 * \note   Copyright (C) 2000-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: CDI.cc 7388 2014-01-22 16:02:07Z kellyt $
//---------------------------------------------------------------------------//

#include "CDI.hh"

namespace rtt_cdi
{
//---------------------------------------------------------------------------//
// Rosseland Spectrum Integrators
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Rosseland spectrum over a frequency range.
 *
 
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Rosseland from x_low to x_high
 *
 */
double CDI::integrateRosselandSpectrum(
    const double lowFreq,
    const double highFreq, 
    const double T)
{
    Require (lowFreq  >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T >= 0.0);

    double planck, rosseland;

    integrate_Rosseland_Planckian_Spectrum(lowFreq, highFreq, T, planck,
                                           rosseland); 

    return rosseland;
}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
// end of CDI_integrateRosselandSpectrum.cc
//---------------------------------------------------------------------------//
