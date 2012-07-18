//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/test/cdi_eospac_test.cc
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:36:36 2001
 * \brief  cdi_eospac test functions.
 * \note   Copyright (C) 2001-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_eospac_test.hh"
#include <iostream>
#include <cmath>

namespace rtt_cdi_eospac_test
{

//---------------------------------------------------------------------------//
// DATA EQUIVALENCE FUNCTIONS USED FOR TESTING
//---------------------------------------------------------------------------//

bool match( const double computedValue,
	    const double referenceValue )
{
    using std::fabs;

    // Start by assuming that the two quantities match exactly.
    bool em = true;
    
    // Compare items up to 10 digits of accuracy.
    
    const double TOL = 1.0e-10;
    
    // Calculate the absolute value of the relative difference between 
    // the computed and reference values.
    
    double reldiff = fabs( ( computedValue - referenceValue )
			   / referenceValue );
    
    // If the comparison fails then change the value of "em" return
    // the result;
    if ( reldiff > TOL )
	em = false;
    
    return em;    
}

//---------------------------------------------------------------------------//

bool match(const std::vector< double >& computedValue, 
	   const std::vector< double >& referenceValue )
{
    using std::fabs;

    // Start by assuming that the two quantities match exactly.
    bool em = true;

    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;

    // Test each item in the list
    double reldiff = 0.0;
    for( size_t i=0; i<computedValue.size(); ++i )
    {
		
	reldiff = fabs( ( computedValue[i] - referenceValue[i] )
			/ referenceValue[i] );
	// If the comparison fails then change the value of "em"
	// and exit the loop.

	if ( reldiff > TOL )
	{
	    em = false;
	    break;
	}
    }
    return em;
}

//---------------------------------------------------------------------------//

bool match(const std::vector< std::vector<double> >& computedValue, 
	   const std::vector< std::vector<double> >& referenceValue )
{
    using std::fabs;

    // Start by assuming that the two quantities match exactly.
    bool em = true;

    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;

    // Test each item in the list
    double reldiff = 0.0;
    for( size_t i=0; i<computedValue.size(); ++i )
    {
	for( size_t j=0; j<computedValue[i].size(); ++j )
	{	    
	    reldiff = fabs( ( computedValue[i][j] - referenceValue[i][j] )
			    / referenceValue[i][j] );

	    // If the comparison fails then change the value of "em"
	    // and exit the loop.
	    if ( reldiff > TOL ) {
		em = false; break; }
	}
    }
    return em;
} 

} // end namespace rtt_cdi_eospac_test

//---------------------------------------------------------------------------//
// end of cdi_eospac_test.cc
//---------------------------------------------------------------------------//


