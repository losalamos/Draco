//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstsign.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:06:56 2004
 * \brief  Test the sign function template.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../sign.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void tstsign( UnitTest & ut )
{    
    if (sign(3.2, 5.6)!=3.2)
    {
	ut.failure("sign: FAILED");
    }
    else
    {
	ut.passes("sign: passed");
    }
    if (sign(4.1, -0.3)!=-4.1)
    {
	ut.failure("sign: FAILED");
    }
    else
    {
	ut.passes("sign: passed");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstsign(ut);
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
//                        end of tstsign.cc
//---------------------------------------------------------------------------//
