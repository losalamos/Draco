//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/Factorial.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide declaration of templatized factorial function.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_factorial_hh
#define sf_factorial_hh

namespace rtt_sf
{

//! \brief factorial
template< typename T >
T factorial( T const k );

//! \brief fraction of factorials, \f$ (k!)/(l!) \f$
template< typename T >
double factorial_fraction( T const k, T const l );

} // end namespace rtt_sf

#endif // sf_Factorial_hh

//---------------------------------------------------------------------------//
//              end of sf/factorial.hh
//---------------------------------------------------------------------------//
