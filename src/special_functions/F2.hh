//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/F2.hh
 * \author Kent Budge
 * \brief  Fermi-Dirac integral of second order.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef sf_F2_hh
#define sf_F2_hh

#include "ds++/config.h"

namespace rtt_sf {
//! Calculate Fermi-Dirac integral of index 2.
DLL_PUBLIC_special_functions double F2(double eta);

} // end namespace rtt_sf

#endif // sf_F2_hh

//---------------------------------------------------------------------------//
// end of sf/F2.hh
//---------------------------------------------------------------------------//
