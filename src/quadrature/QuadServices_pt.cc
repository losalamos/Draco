//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices_pt.cc
 * \author Kelly Thompson
 * \date   Mon Nov  8 14:23:03 2004
 * \brief  Member definitions of class QuadServices
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "QuadServices.i.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
// Make kronecker delta is valid for double, int, unsigned, float

template 
double QuadServices::kronecker_delta( double const test_value, 
				      double const offset ) const;
template 
int QuadServices::kronecker_delta( int const test_value, 
				   int const offset ) const;
template 
unsigned QuadServices::kronecker_delta( unsigned const test_value, 
					unsigned const offset ) const;
template 
float QuadServices::kronecker_delta( float const test_value, 
				     float const offset ) const;

//---------------------------------------------------------------------------//
// Make factorial valid for int, unsigned    

template
unsigned QuadServices::factorial( unsigned const k ) const;

template
int QuadServices::factorial( int const k ) const;

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices_pt.cc
//---------------------------------------------------------------------------//
