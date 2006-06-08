//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/Ylm.hh
 * \author Kent Budge
 * \date   Tue Jul  6 10:03:25 MDT 2004
 * \brief  Declare the Ylm function template.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef special_functions_Ylm_hh
#define special_functions_Ylm_hh

namespace rtt_sf
{

//! Compute the spherical harmonic \f$ Y_{lm}(\theta,\phi) \f$ 
double Ylm(unsigned l, int m, double theta, double phi);

} // end namespace rtt_sf

#endif // special_functions_Ylm_hh

//---------------------------------------------------------------------------//
//              end of utils/Ylm.hh
//---------------------------------------------------------------------------//
