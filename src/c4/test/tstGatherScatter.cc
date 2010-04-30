//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstGatherScatter.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::gather and c4::scatter functions
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "../Release.hh"
#include "../ParallelUnitTest.hh"
#include "../C4_Functions.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstDeterminateGatherScatter(UnitTest &ut)
{
    unsigned pid = node();
    unsigned const number_of_processors = nodes();
    vector<unsigned> pids(number_of_processors);
    gather(&pid, &pids[0], 1);

    if (pid == 0)
    {
        for (unsigned i=1; i<number_of_processors; ++i)
        {
            pids[i] += pids[i-1];
        }
    }

    unsigned base;
    scatter(&pids[0], &base, 1);

    if (base == pid*(pid+1)/2)
    {
        ut.passes("correct base summation in gather/scatter");
    }
    else
    {
        ut.failure("NOT correct base summation in gather/scatter");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstDeterminateGatherScatter(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstGatherScatter, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstGatherScatter, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstGatherScatter.cc
//---------------------------------------------------------------------------//
