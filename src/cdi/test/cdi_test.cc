//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/cdi_test.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 10:51:39 2001
 * \brief  CDI Test help functions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_test.hh"
#include <iostream>
#include <cmath>

namespace rtt_cdi_test
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
// CHECK COMPUTED VERSUS EXPECTED VALUES
//---------------------------------------------------------------------------//

bool match( double computedValue, double referenceValue )
{
    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;
	    
    // Calculate the absolute value of the relative difference between 
    // the computed and reference values.
    double reldiff = std::fabs( ( computedValue - referenceValue )
				/ referenceValue );
    
    // If the comparison fails then return "false" to indicate that
    // the test failed.
    if ( reldiff > TOL ) return false;
    
    return true;        
} 

bool match( const std::vector< double > &computedValue, 
	    const std::vector< double > &referenceValue ) 
{
    // If the vector sizes don't match then die
    if ( computedValue.size() != referenceValue.size() )
	return false;
	    
    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;
	    
    // Test each item in the list
    double reldiff = 0.0;
    for ( int i=0; i<computedValue.size(); ++i )
    {
	reldiff = std::fabs( ( computedValue[i] - referenceValue[i] )
			     / referenceValue[i] );
	// If the comparison fails then return "false" to indicate 
	// that the test failed.
	if ( reldiff > TOL ) return false;
    }
    return true;
}

} // end namespace rtt_cdi_test

//---------------------------------------------------------------------------//
//                              end of cdi_test.cc
//---------------------------------------------------------------------------//
