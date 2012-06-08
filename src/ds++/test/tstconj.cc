//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstconj.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <complex>

#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../conj.hh"
#include "../square.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstconj(UnitTest &ut)
{
    if (soft_equiv(conj(3.5), 3.5))
    {
	ut.passes("conj(double) is correct");
    }
    else
    {
	ut.failure("conj(double) is NOT correct");
    }

    complex<double> c(2.7, -1.4);
    if (soft_equiv((rtt_dsxx::conj(c)*c).real(), square(abs(c))))
    {
	ut.passes("conj(complex) is correct");
    }
    else
    {
	ut.failure("conj(complex) is NOT correct");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc,argv,release );
	tstconj(ut);
    }
    catch (exception &err)
    {
	cout << "ERROR: While testing tstconj, " << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
	cout << "ERROR: While testing tstconj, " 
             << "An unknown exception was thrown." << endl;
	return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstconj.cc
//---------------------------------------------------------------------------//
