//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstField_Traits.cc
 * \author Kent Budge
 * \date   Tue Aug 26 12:18:55 2008
 * \brief  
 * \note   Copyright (C) 2006-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../Field_Traits.hh"
#include "../value.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        if (Field_Traits<complex<double> >::zero() == 0.0)
            ut.passes("complex zero good");
        else
            ut.failure("complex zero NOT good");
        if (Field_Traits<complex<double> >::one() == 1.0)
            ut.passes("complex zero good");
        else
            ut.failure("complex zero NOT good");
        double const x = 3.7;
        if (value(x) == 3.7)
            ut.passes("complex zero good");
        else
            ut.failure("complex zero NOT good");
        
        if (Field_Traits<double const>::zero() == 0.0)
            ut.passes("double zero good");
        else
            ut.failure("double zero NOT good");
        if (Field_Traits<double const>::one() == 1.0)
            ut.passes("double zero good");
        else
            ut.failure("double zero NOT good");
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstField_Traits.cc
//---------------------------------------------------------------------------//
