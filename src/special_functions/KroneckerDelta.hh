//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/KroneckerDelta.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide declaration of templatized KroneckerDelta function.
 * \note   Copyright 2006 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_KroneckerDelta_hh
#define sf_KroneckerDelta_hh

namespace rtt_sf
{

//! \brief kronecker_delta
template< typename T >
T kronecker_delta( T const test_value, T const offset );

} // end namespace rtt_sf

#endif // sf_KroneckerDelta_hh

//---------------------------------------------------------------------------//
//              end of sf/KroneckerDelta.hh
//---------------------------------------------------------------------------//
