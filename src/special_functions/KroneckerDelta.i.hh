//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/KroneckerDelta.i.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide implementation of templatized delta function.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_KroneckerDelta_i_hh
#define sf_KroneckerDelta_i_hh

#include "KroneckerDelta.hh"

namespace rtt_sf {

//---------------------------------------------------------------------------//
/*! 
 * \brief kronecker_delta
 * 
 * Return 1 if test_value == offset, otherwise return 0;
 * 
 * \param test_value
 * \param offset
 * \return 1 if test_value == offset, otherwise return 0;
 */
template <typename T>
unsigned int kronecker_delta(T const test_value, T const offset) {
  return (test_value == offset) ? 1 : 0;
}

} // end namespace rtt_sf

#endif // sf_KroneckerDelta_i_hh

//---------------------------------------------------------------------------//
//              end of sf/KroneckerDelta.i.hh
//---------------------------------------------------------------------------//
