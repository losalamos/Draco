//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstdbc.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Copyright (C) 2003-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../isFinite.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void isFinite_test( UnitTest & ut )
{
    using rtt_dsxx::isFinite;
    
    double x(15.0);
    if( isFinite(x) )
    {
	ut.passes("Correctly found 15. to be finite.");
    }
    else
    {
	ut.failure("Failed to find 15. to be finite.");
    }
    
    if( !isNaN(x) )
    {
	ut.passes("Correctly found 15. to be not NaN.");
    }
    else
    {
	ut.failure("Incorrectly found 15. to be finite.");
    }
   
    if( !isInfinity(x) )
    {
	ut.passes("Correctly found 15. to be not infinity.");
    }
    else
    {
	ut.failure("Incorrectly found 15. to be infinity.");
    }
        
    double y(x/(15.0-x));
    if( isInfinity(y) )
    {
	ut.passes("Correctly found 15./0. to be infinity.");
    }
    else
    {
	ut.failure("Incorrectly found 15./0. to be finite.");
    }

    if( !isNaN(y) )
    {
	ut.passes("Correctly found 15./0. to be not NaN.");
    }
    else
    {
	ut.failure("Incorrectly found 15./0. to be NaN.");
    }
    
    if( !isFinite(y) )
    {
	ut.passes("Correctly found 15./0. to be not finite.");
    }
    else
    {
	ut.failure("Incorrectly found 15./0. to be finite.");
    }
       
    double z(0.0*y);
    if(!isInfinity(z) )
    {
	ut.passes("Correctly found 0.*infinity to be not infinite.");
    }
    else
    {
	ut.failure("Incorrectly found 0*infinity to be infinite.");
    }

    if( isNaN(z) )
    {
	ut.passes("Correctly found 0.*infinity to be NaN.");
    }
    else
    {
	ut.failure("Failed to find 0.*infinity to be NaN.");
    }
    
    if( !isFinite(z) )
    {
	ut.passes("Correctly found 0.*infinity to be not finite.");
    }
    else
    {
	ut.failure("Incorrectly found 0.*infinity to be finite.");
    }
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut( argc, argv, release );
    try
    {
	isFinite_test(ut);
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstdbc.cc
//---------------------------------------------------------------------------//
