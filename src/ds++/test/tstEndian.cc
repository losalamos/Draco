//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstEndian.cc
 * \author Mike Buksas
 * \date   Tue Oct 23 16:20:59 2007
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>

#include "ds++/Assert.hh"
#include "ds++/ScalarUnitTest.hh"
#include "../Release.hh"

#include "../Endian.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_char_data(ScalarUnitTest& ut)
{

    unsigned char data[] = {'a', 'b', 'c'};
    unsigned int length = sizeof(data)/sizeof(unsigned char);

    char_byte_swap(data, length);

    if ((data[0] != 'c') || (data[1] != 'b') || (data[2] != 'a'))
        ut.failure("char_byte_swap function failed");
    
}

void test_integer(ScalarUnitTest& ut)
{

    // Integer. This value overflows unsigned ints.
    int moo = 0xDEADBEEF;

    byte_swap(moo);

    if (moo != 0xEFBEADDE)
        ut.failure("byte_swap failed for for integer type");

    // Unsigned integer
    unsigned int u_moo = 0xDEADBEEF;

    byte_swap(u_moo);

    if (u_moo != 0xEFBEADDE)
        ut.failure("byte_swap failed for for unsigned integer type");

}


void test_idempotence(ScalarUnitTest& ut)
{
    
    /* This test demonstrates that two applications of byte-swap in succession
     * return the original value.
     *
     * To do this, we sweep over a lot of double values. We use a non-integral
     * multiplier for successive values to avoid small subsets of the
     * available patterns of bits. E.g. multiples of 2.
     */ 

    for (double value = 1.0;
         value < std::numeric_limits<double>::max();
         value *= 3.4)
    {

        double local = value;
        byte_swap(local);   
        byte_swap(local);

        // These numbers should be identical, so I'm testing for equality.
        if (local != value)
            ut.failure("byte_swap failed to reproduce original number");


        double neg_local = -value;
        byte_swap(neg_local);
        byte_swap(neg_local);

        if (neg_local != -value)
            ut.failure("byte_swap failed to reproduce original number");

    }

}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, release);
    try
    {
        test_char_data(ut);
        test_integer(ut);
        test_idempotence(ut);
        ut.passes("Just Because.");
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstEndian, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstEndian, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstEndian.cc
//---------------------------------------------------------------------------//
