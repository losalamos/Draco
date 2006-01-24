//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/test/test_sf.cc
 * \author Kelly Thompson
 * \date   Tue Sep 27 12:49:39 2005
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "ds++/Assert.hh"
#include "../Factorial.hh"
#include "../KroneckerDelta.hh"
#include "../Release.hh"
#include "sf_test.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstKdelta()
{
    using rtt_sf::kronecker_delta;
    if( kronecker_delta( 0, 0 ) == 1 )
    {
	PASSMSG("Found kronecker_delta(0,0) == 1, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta(0,0) != 1, kronecker_delta is not working.");
    }
    if( kronecker_delta( 0, 1 ) == 0 )
    {
	PASSMSG("Found kronecker_delta(0,1) == 0, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta(0,1) != 0, kronecker_delta is not working.");
    }
    if( kronecker_delta( 1, 1 ) == 1 )
    {
	PASSMSG("Found kronecker_delta(1,1) == 1, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta(1,1) != 1, kronecker_delta is not working.");
    }
    if( kronecker_delta( 1, 0 ) == 0 )
    {
	PASSMSG("Found kronecker_delta(1,0) == 0, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta(1,0) != 0, kronecker_delta is not working.");
    }
    if( kronecker_delta( -1, 0 ) == 0 )
    {
	PASSMSG("Found kronecker_delta(-1,0) == 0, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta(-1,0) != 0, kronecker_delta is not working.");
    }
    if( kronecker_delta( -1, -1 ) == 1 )
    {
	PASSMSG("Found kronecker_delta(-1,-1) == 1, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta(-1,-1) != 1, kronecker_delta is not working.");
    }

    unsigned uZero(0);
    unsigned uOne(1);
    long lZero(0);
    long lOne(1);

    if( kronecker_delta( uOne, uZero ) == uZero )
    {
	PASSMSG("Found kronecker_delta<unsigned>(uOne,uZero) == uZero, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta<unsigned>(uOne,uZero) != uZero, kronecker_delta is not working.");
    }

    if( kronecker_delta( lOne, lZero ) == lZero )
    {
	PASSMSG("Found kronecker_delta<long>(uOne,uZero) == uZero, kronecker_delta is working.");
    }
    else
    {
	PASSMSG("Found kronecker_delta<long>(uOne,uZero) != uZero, kronecker_delta is not working.");
    }

    return;
}

//---------------------------------------------------------------------------//

void tstFactorial()
{
    using rtt_sf::factorial;

    // Test factorial

    if( factorial(0) == 1 )
    {
	PASSMSG("Found factorial(0) == 1, factorial is working.");
    }
    else
    {
	PASSMSG("Found factorial(0) != 1, factorial is not working.");
    }
    if( factorial(1) == 1 )
    {
	PASSMSG("Found factorial(1) == 1, factorial is working.");
    }
    else
    {
	PASSMSG("Found factorial(1) != 1, factorial is not working.");
    }
    if( factorial(2) == 2 )
    {
	PASSMSG("Found factorial(2) == 2, factorial is working.");
    }
    else
    {
	PASSMSG("Found factorial(2) != 2, factorial is not working.");
    }
    if( factorial(3) == 6 )
    {
	PASSMSG("Found factorial(3) == 6, factorial is working.");
    }
    else
    {
	PASSMSG("Found factorial(3) != 6, factorial is not working.");
    }
    if( factorial(-3) == 1 )
    {
	PASSMSG("Found factorial(-3) == 1, factorial is working.");
    }
    else
    {
	PASSMSG("Found factorial(-3) != 1, factorial is not working.");
    }

    unsigned uOne(1);
    long     lOne(1);
    
    if( factorial(uOne) == uOne )
    {
	PASSMSG("Found factorial<unsigned>(1) == unsigned(1), factorial is working.");
    }
    else
    {
	PASSMSG("Found factorial<unsigned>(1) != unsigned(1), factorial is not working.");
    }
    if( factorial(lOne) == lOne )
    {
	PASSMSG("Found factorial<long>(1) == long(1), factorial is working.");
    }
    else
    {
	PASSMSG("Found factorial<long>(1) != long(1), factorial is not working.");
    }   
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;
    using std::string;
    
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (string(argv[arg]) == "--version")
        {
            cout << argv[0] << ": version " << rtt_sf::release() << endl;
            return 0;
        }

    try
    {
        tstFactorial();
        tstKdelta();
    }
    catch (std::exception &err)
    {
        cout << "ERROR: While testing test_sf," 
                  << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing test_sf, " 
                  << "An unknown exception was thrown." << endl;
        return 1;
    }

    // status of test
    cout <<   "\n*********************************************\n";
    if (rtt_sf_test::passed) 
        cout << "**** test_sf Test: PASSED";
    cout <<   "\n*********************************************\n\n"
         << "Done testing test_sf." << endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of test_sf.cc
//---------------------------------------------------------------------------//
