//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI_integrate_Rosseland_Planckian_Spectrum.cc
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
 * \brief Integrate the Planckian and Rosseland spectrum over a frequency
 *        range.
 * \param T the temperature in keV (must be greater than 0.0)
 * \return void the integrated normalized Planckian and Rosseland from x_low
 *        to x_high are passed by reference in the function call
 *
 * \f[
 * planck(T) = \int_{\nu_1}^{\nu_2}{B(\nu,T)d\nu}
 * rosseland(T) = \int_{\nu_1}^{\nu_2}{\frac{\partial B(\nu,T)}{\partial T}d\nu}
 * \f]
 */
void CDI::integrate_Rosseland_Planckian_Spectrum(double lowFreq,
						 double highFreq,
						 double const T,
						 double &planck, 
						 double &rosseland)
{
    Require (lowFreq >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T >= 0.0);

    if (T==0.0)
    {
        planck = 0.0;
        rosseland = 0.0;
        return;
    }

    double planck_high, rosseland_high;
    double planck_low,  rosseland_low;

    // Sale the frequencies by temperature
    lowFreq /= T;
    highFreq /= T;

    double const exp_lowFreq  = std::exp(-lowFreq);
    double const exp_highFreq = std::exp(-highFreq);

    integrate_planck_rosseland(lowFreq,  exp_lowFreq,  planck_low,
                               rosseland_low);
    integrate_planck_rosseland(highFreq, exp_highFreq, planck_high,
                               rosseland_high);

    planck    = planck_high    - planck_low;
    rosseland = rosseland_high - rosseland_low;

    return;
}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
// end of CDI_integrate_Rosseland_Planckian_Spectrum.cc
//---------------------------------------------------------------------------//
