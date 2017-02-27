//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/KroneckerDelta.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide declaration of templatized KroneckerDelta function.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_KroneckerDelta_hh
#define sf_KroneckerDelta_hh

#include "ds++/config.h"

namespace rtt_sf {

//! \brief kronecker_delta
template <typename T>
DLL_PUBLIC_special_functions unsigned int kronecker_delta(T const test_value,
                                                          T const offset);

} // end namespace rtt_sf

#endif // sf_KroneckerDelta_hh

//---------------------------------------------------------------------------//
// end of sf/KroneckerDelta.hh
//---------------------------------------------------------------------------//
