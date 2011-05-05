//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   roots/test/tstquadratic.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "ds++/Release.hh"
#include "../quadratic.hh"

using namespace std;
using namespace rtt_roots;
using namespace rtt_dsxx;


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//



//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );

        double r1, r2;
        quadratic(2.0, -2.0, -12.0, r1, r2);

        if (r1>r2)
        {
            swap(r1, r2);
        }
        
	if (soft_equiv(r1, -2.0))
	{
	    ut.passes("quadratic solve returned correct first root");
	}
	else
	{
	    ut.failure("quadratic solve returned INCORRECT first root");
	}
        
	if (soft_equiv(r2, 3.0))
	{
	    ut.passes("quadratic solve returned correct second root");
	}
	else
	{
	    ut.failure("quadratic solve returned INCORRECT second root");
	}
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown." << endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of testquadratic.cc
//---------------------------------------------------------------------------//
