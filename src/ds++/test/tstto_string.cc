//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstto_string.cc
 * \author Kent Budge
 * \date   Fri Jul 25 08:49:48 2008
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <ds++/config.h> // Must be first!

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "../Assert.hh"
#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../to_string.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstto_string( UnitTest &ut)
{
    string const pi = to_string(M_PI, 20);

    if (soft_equiv(M_PI, atof(pi.c_str())))
    {
        ut.passes("pi correctly written/read");
    }
    else
    {
        ut.failure("pi NOT correctly written/read");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstto_string(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstto_string, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstto_string, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstto_string.cc
//---------------------------------------------------------------------------//
