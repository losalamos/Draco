//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   special_functions/F12inv.hh
 * \author Kent Budge
 * \date   Tue Sep 21 09:20:10 2004
 * \brief  Inverse Fermi-Dirac integral of 1/2 order.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef sf_F12inv_hh
#define sf_F12inv_hh

#include "ds++/config.h"

namespace rtt_sf {
//! Compute the inverse Fermi-Dirac function of index 1/2.
DLL_PUBLIC_special_functions double F12inv(double f);

//! Compute the inverse Fermi-Dirac function of index 1/2 and its derivative.
DLL_PUBLIC_special_functions void F12inv(double f, double &eta, double &deta);

} // end namespace rtt_sf

#endif // sf_F12inv_hh

//----------------------------------------------------------------------------//
// end of sf/F12inv.hh
//----------------------------------------------------------------------------//
