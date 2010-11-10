//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstcube.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../cube.hh"

using namespace std;
using namespace rtt_dsxx;
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
	if (soft_equiv(cube(2.0), 8.0))
	{
	    ut.passes("square function returned correct double");
	}
	else
	{
	    ut.failure("square function did NOT return correct double.");
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
//                        end of testcube.cc
//---------------------------------------------------------------------------//
