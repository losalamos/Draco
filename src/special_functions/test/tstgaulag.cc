//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/test/tstgaulag.cc
 * \author Kent Budge
 * \date   Tue Sep 27 12:49:39 2005
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id: tstgaulag.cc 5830 2011-05-05 19:43:43Z kellyt $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"
#include "ds++/square.hh"
#include "ds++/ScalarUnitTest.hh"
#include "../gaulag.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_special_functions;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstgaulag(UnitTest &ut)
{
    vector<double> x, w;
    gaulag(x, w, 0.0, 3);
    double sum = 0.0;
    for (unsigned i=0; i<3; ++i)
    {
        sum += x[i]*square(square(x[i]))*w[i];
    }
    if (!soft_equiv(sum, 120.0))
    {
        ut.failure("gaulag NOT accurate");
    }
    else
    {
        ut.passes("gaulag accurate");
    }        
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut(argc, argv, release);
        tstgaulag(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstgaulag, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstgaulag, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstgaulag.cc
//---------------------------------------------------------------------------//
