//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/test/cdi_eospac_test.cc
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:36:36 2001
 * \brief  cdi_eospac test functions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_eospac_test.hh"
#include <iostream>
#include <cmath>

namespace rtt_cdi_eospac_test
{

//===========================================================================//
// PASS/FAILURE
//===========================================================================//

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
	      << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool pass_msg(const std::string &passmsg)
{
    std::cout << "Test: passed" << std::endl;
    std::cout << "     " << passmsg << std::endl;
    return true;
}

//---------------------------------------------------------------------------//

bool fail_msg(const std::string &failmsg)
{
    std::cout << "Test: failed" << std::endl;
    std::cout << "     " << failmsg << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//
// BOOLEAN PASS FLAG
//---------------------------------------------------------------------------//

bool passed = true;

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
    for ( int i=0; i<computedValue.size(); ++i )
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
    for ( int i=0; i<computedValue.size(); ++i )
    {
	for ( int j=0; j<computedValue[i].size(); ++j )
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


