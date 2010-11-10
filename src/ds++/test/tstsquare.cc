//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSQR.cc
 * \author Kent Budge
 * \date   Tue Jul  6 08:54:38 2004
 * \brief  Test the SQR function.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../square.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//


int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	if (soft_equiv(square(3.0), 9.0))
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
//                        end of tstSQR.cc
//---------------------------------------------------------------------------//
