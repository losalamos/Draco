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

#include "QuadServices.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
// Make kronecker delta is valid for double, int, unsigned, float

template 
double kronecker_delta( double const test_value, 
			double const offset ) ;
template 
int kronecker_delta( int const test_value, 
		     int const offset ) ;
template 
unsigned kronecker_delta( unsigned const test_value, 
			  unsigned const offset ) ;
template 
float kronecker_delta( float const test_value, 
		       float const offset ) ;

//---------------------------------------------------------------------------//
// Make factorial valid for int, unsigned    

template
unsigned factorial( unsigned const k ) ;

template
int factorial( int const k ) ;

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices_pt.cc
//---------------------------------------------------------------------------//
