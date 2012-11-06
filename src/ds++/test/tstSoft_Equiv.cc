//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSoft_Equiv.cc
 * \author Thomas M. Evans
 * \date   Wed Nov  7 15:55:54 2001
 * \brief  Soft_Equiv header testing utilities.
 * \note   Copyright (C) 2001-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/config.h"
#include "../Release.hh"
#include "../Soft_Equivalence.hh"
#include "../Assert.hh"
#include "../ScalarUnitTest.hh"
#include <vector>
#include <deque>
#ifdef HAS_CXX11_ARRAY
#include <array>
#endif

#define ITFAILS    ut.failure( __LINE__, __FILE__ )
#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)

using namespace std;
using rtt_dsxx::soft_equiv;
using rtt_dsxx::soft_equiv_deep;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_soft_equiv_scalar(rtt_dsxx::ScalarUnitTest & ut)
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
    return;
}

//---------------------------------------------------------------------------//

void test_soft_equiv_container(rtt_dsxx::ScalarUnitTest & ut)
{
    vector<double> values(3, 0.0);
    values[0] = 0.3247333291470;
    values[1] = 0.3224333221471;
    values[2] = 0.3324333522912;
    
    vector<double> const reference( values );

    if (soft_equiv(values.begin(), values.end(),
		   reference.begin(), reference.end()))
	PASSMSG("Passed vector equivalence test.");
    else
	ITFAILS;

    // modify one value (delta < tolerance )
    values[1] += 1.0e-13;
    if (!soft_equiv(values.begin(), values.end(),
		    reference.begin(), reference.end(), 1.e-13))
	PASSMSG("Passed vector equivalence precision test.");
    else
	ITFAILS;

    // Tests that compare 1D vector data to 1D array data.
    double v[3];
    v[0] = reference[0];
    v[1] = reference[1];
    v[2] = reference[2];

    if (soft_equiv(&v[0], &v[3],
		   reference.begin(), reference.end()))    
	PASSMSG("Passed vector-pointer equivalence test.");
    else
	ITFAILS;

    if (!soft_equiv(reference.begin(), reference.end(), &v[0], &v[3]))
	ITFAILS;

    // modify one value (delta < tolerance )
    v[1] += 1.0e-13;
    if (!soft_equiv(&v[0], v+3,
		    reference.begin(), reference.end(), 1.e-13))
	PASSMSG("Passed vector-pointer equivalence precision test.");
    else
	ITFAILS;

#ifdef HAS_CXX11_ARRAY
#ifdef HAS_CXX11_INITIALIZER_LISTS
    // C++ std::array containers
    std::array<double,3> cppa_vals{
        { 0.3247333291470, 0.3224333221471, 0.3324333522912 } };
    if (soft_equiv(cppa_vals.begin(), cppa_vals.end(),
		   reference.begin(), reference.end()))
	PASSMSG("Passed std::array<int,3> equivalence test.");
    else
	ITFAILS;
#endif
#endif

    // Try with a std::deque
    deque<double> d;
    d.push_back(reference[0]);
    d.push_back(reference[1]);
    d.push_back(reference[2]);
    if (soft_equiv(d.begin(), d.end(),
                   reference.begin(), reference.end()))
	PASSMSG("Passed deque<T> equivalence test.");
    else
	ITFAILS;
    
    return;
}

//---------------------------------------------------------------------------//

#ifdef HAS_CXX11_ARRAY
#ifdef HAS_CXX11_INITIALIZER_LISTS
void test_soft_equiv_deep_container(rtt_dsxx::ScalarUnitTest & ut)
{
    
    vector<vector<double> > values = {
        { 0.3247333291470, 0.3224333221471, 0.3324333522912 },
        { 0.3247333292470, 0.3224333222471, 0.3324333523912 },
        { 0.3247333293470, 0.3224333223471, 0.3324333524912 }
    };
    vector<vector<double> > const reference = values; 
        
    if (soft_equiv_deep<2>().equiv(values.begin(), values.end(),
                                   reference.begin(), reference.end()))
        PASSMSG("Passed vector<vector<double>> equivalence test.");
    else
        ITFAILS;

    // Soft_Equiv should still pass
    values[0][1]+= 1.0e-13;
    if (!soft_equiv_deep<2>().equiv(values.begin(), values.end(),
                                    reference.begin(), reference.end(), 1.e-13))
        PASSMSG("Passed vector<vector<double>> equivalence precision test.");
    else
        ITFAILS;

    // Compare C++ array to vector<vector<double>> data.
    // This cannot work because the C++ array is fundamentally a 1-D container.
        
    // double v[3];
    // v[0] = 0.3247333291470;
    // v[1] = 0.3224333221471;
    // v[2] = 0.3324333522912;
    // if (soft_equiv(&v[0], &v[3],
    //     	   reference.begin(), reference.end()))    
    //     PASSMSG("Passed vector-pointer equivalence test.");
    // else
    //     ITFAILS;
       
    // if (!soft_equiv(reference.begin(), reference.end(), &v[0], &v[3]))
    //     ITFAILS;

    // Test 3-D array
    vector<vector<vector<double> > > const ref = {
        {
            { 0.1, 0.2 },
            { 0.3, 0.4 }
        },
        {
            { 1.1, 1.2 },
            { 1.3, 1.4 }
        },
        {
            { 2.1, 2.2 },
            { 2.3, 2.4 }
        }
    };
    vector<vector<vector<double> > > val = ref;
    
    if (soft_equiv_deep<3>().equiv(val.begin(), val.end(),
                                   ref.begin(), ref.end()))
        PASSMSG("Passed vector<vector<vector<double>>> equivalence test.");
    else
        ITFAILS;

    return;
}
#endif
#endif

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        // >>> UNIT TESTS
	test_soft_equiv_scalar(ut);
	test_soft_equiv_container(ut);
#ifdef HAS_CXX11_ARRAY
#ifdef HAS_CXX11_INITIALIZER_LISTS
	test_soft_equiv_deep_container(ut);
#endif
#endif
    }
    catch (rtt_dsxx::assertion &excpt)
    {
        std::cout << "ERROR: While testing tstSoft_Equiv, "
                  << excpt.what() << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSoft_Equiv, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
// end of tstSoft_Equiv.cc
//---------------------------------------------------------------------------//
