//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI_integratePlanckSpectrum.cc
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
/*!
 * \brief Integrate the Planckian spectrum over a frequency range.
 *
 * \param lowFreq lower frequency bound in keV
 * \param highFreq higher frequency bound in keV
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian from x_low to x_high
 *
 */
double CDI::integratePlanckSpectrum(double lowFreq, 
				    double highFreq, 
				    const double T) 
{
    Require (lowFreq  >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T        >= 0.0);

    // return 0 if temperature is a hard zero
    if (T == 0.0) return 0.0;

    // Sale the frequencies by temperature
    lowFreq /= T;
    highFreq /= T;

    double integral =
        integrate_planck(highFreq) - integrate_planck(lowFreq);
    

    Ensure ( integral >= 0.0 );
    Ensure ( integral <= 1.0 );

    return integral;
}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
// end of CDI_integratePlanckSpectrum.cc
//---------------------------------------------------------------------------//
