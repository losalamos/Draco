//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/KroneckerDelta_pt.cc
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide explicit instantiations of templatized KroneckerDelta function.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "KroneckerDelta.i.hh"

namespace rtt_sf {

//---------------------------------------------------------------------------//
// Make kronecker delta valid only for double, int, unsigned, and float.

template DLL_PUBLIC_special_functions unsigned int
kronecker_delta(double const test_value, double const offset);
template DLL_PUBLIC_special_functions unsigned int
kronecker_delta(int const test_value, int const offset);
template DLL_PUBLIC_special_functions unsigned int
kronecker_delta(long const test_value, long const offset);
template DLL_PUBLIC_special_functions unsigned int
kronecker_delta(unsigned const test_value, unsigned const offset);
template DLL_PUBLIC_special_functions unsigned int
kronecker_delta(float const test_value, float const offset);

} // end namespace rtt_sf

//---------------------------------------------------------------------------//
// end of sf/KroneckerDelta_pt.cc
//---------------------------------------------------------------------------//
