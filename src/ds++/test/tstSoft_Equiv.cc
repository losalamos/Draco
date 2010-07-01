//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSoft_Equiv.cc
 * \author Thomas M. Evans
 * \date   Wed Nov  7 15:55:54 2001
 * \brief  Soft_Equiv header testing utilities.
 * \note   Copyright (c) 2001-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds_test.hh"
#include "../Release.hh"
#include "../Soft_Equivalence.hh"
#include "../Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <list>

using namespace std;

using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_soft_equiv_scalar()
{
    // ensure that we can not use integer tolerance.
    {
        int x = 31415;
        int y = 31416;
        int tol = 1l;

        try
        {
            /* bool result =  */ soft_equiv(x,y,tol);
            throw "Bogus!";            
        }
        catch( rtt_dsxx::assertion const & /* error */ )
        {
            PASSMSG("Successfully prevented use of soft_equiv(int,int,int).");
        }
        catch( ... )
        {
            FAILMSG("We should never get here.");
        }
    }
    
    // test with doubles
    {
        double x = 0.9876543212345678;
        double y = 0.9876543212345678;
        
        if (!soft_equiv(x, y, 1.e-16)) ITFAILS;
        if (!soft_equiv(x, y))         ITFAILS;
        
        double z = 0.9876543212345679;
        
        if (soft_equiv(x, z, 1.e-16)) ITFAILS;
        
        double a = 0.987654321234;
        
        if (!soft_equiv(x, a)) ITFAILS;
        
        a = 0.987654321233;
        
        if (soft_equiv(x, a)) ITFAILS;      

        // checks for the new "reference=zero" coding 4aug00
        double zero = 0.0;
        if ( soft_equiv( 1.0e-10, zero)) ITFAILS;
        if ( soft_equiv(-1.0e-10, zero)) ITFAILS;
        if (!soft_equiv(-1.0e-35, zero)) ITFAILS;
        if (!soft_equiv( 1.0e-35, zero)) ITFAILS;
    }
}

//---------------------------------------------------------------------------//

void test_soft_equiv_container()
{
    vector<double> values(3, 0.0);
    values[0] = 0.3247333291470;
    values[1] = 0.3224333221471;
    values[2] = 0.3324333522912;
    
    vector<double> reference(3, 0.0);
    reference[0] = 0.3247333291470;
    reference[1] = 0.3224333221471;
    reference[2] = 0.3324333522912;

    if (soft_equiv(values.begin(), values.end(),
		   reference.begin(), reference.end()))
    {
	PASSMSG("Passed vector equivalence test.");
    }
    else
    {
	ITFAILS;
    }


    values[1] = 0.3224333221472;
    if (!soft_equiv(values.begin(), values.end(),
		    reference.begin(), reference.end(), 1.e-13))
    {
	PASSMSG("Passed vector equivalence precision test.");
    }
    else
    {
	ITFAILS;
    }

    double v[3];
    v[0] = 0.3247333291470;
    v[1] = 0.3224333221471;
    v[2] = 0.3324333522912;

    if (soft_equiv(&v[0], &v[3],
		   reference.begin(), reference.end()))    
    {
	PASSMSG("Passed vector-pointer equivalence test.");
    }
    else
    {
	ITFAILS;
    }

    if (!soft_equiv(reference.begin(), reference.end(), &v[0], &v[3]))
	ITFAILS;

    v[1] = 0.3224333221472;
    if (!soft_equiv(&v[0], v+3,
		    reference.begin(), reference.end(), 1.e-13))
    {
	PASSMSG("Passed vector-pointer equivalence precision test.");
    }
    else
    {
	ITFAILS;
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	test_soft_equiv_scalar();
	test_soft_equiv_container();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstSoft_Equiv, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_ds_test::passed) 
    {
        cout << "**** tstSoft_Equiv Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstSoft_Equiv." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstSoft_Equiv.cc
//---------------------------------------------------------------------------//
