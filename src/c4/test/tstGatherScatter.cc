//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstGatherScatter.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::gather and c4::scatter functions
 * \note   Copyright (C) 2006-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "../ParallelUnitTest.hh"
#include "../global.hh"
#include "../gatherv.hh"
#include "../scatterv.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <vector>
#include <cmath>

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
void tstIndeterminateGatherScatterv(UnitTest &ut)
{
    unsigned pid = node();
    unsigned const number_of_processors = nodes();
    vector<unsigned> send(pid, pid);
    vector<vector<unsigned> > receive;
    indeterminate_gatherv(send, receive);

    ut.passes("No exception thrown");

    if (pid==0)
    {
        if (receive.size() == number_of_processors)
        {
            ut.passes("correct number of processors in gatherv");
        }
        else
        {
            ut.failure("NOT correct number of processors in gatherv");
        }
        for (unsigned p=0; p<number_of_processors; ++p)
        {
            if (receive[p].size() != p)
            {
                ut.failure("NOT correct number of elements in gatherv");
                for (unsigned i=0; i<p; ++i)
                {
                    if (receive[p][i] != p)
                    {
                        ut.failure("NOT correct values in gatherv");
                    }
                }
            }
            // Prepare for next test
            receive[p].resize(0);
            receive[p].resize(2*p, 3*p);
        }
    }

    indeterminate_scatterv(receive, send);

    if (send.size() == 2*pid)
    {
        ut.passes("correct number of processors in scatterv");
    }
    else
    {
        ut.failure("NOT correct number of processors in scatterv");
    }
    for (unsigned i=0; i<2*pid; ++i)
    {
        if (send[i] != 3*pid)
        {
            ut.failure("NOT correct values in scatterv");
        }
    }
}

//---------------------------------------------------------------------------//
void tstDeterminateGatherScatterv(UnitTest &ut)
{
    unsigned pid = node();
    unsigned const number_of_processors = nodes();
    vector<unsigned> send(pid, pid);
    vector<vector<unsigned> > receive(number_of_processors);
    for (unsigned p=0; p<number_of_processors; ++p)
    {
        receive[p].resize(p, p);
    }
    determinate_gatherv(send, receive);

    ut.passes("No exception thrown");

    if (pid==0)
    {
        if (receive.size() == number_of_processors)
        {
            ut.passes("correct number of processors in gatherv");
        }
        else
        {
            ut.failure("NOT correct number of processors in gatherv");
        }
        for (unsigned p=0; p<number_of_processors; ++p)
        {
            if (receive[p].size() != p)
            {
                ut.failure("NOT correct number of elements in gatherv");
                for (unsigned i=0; i<p; ++i)
                {
                    if (receive[p][i] != p)
                    {
                        ut.failure("NOT correct values in gatherv");
                    }
                }
            }
            // Prepare for next test
            receive[p].resize(0);
            receive[p].resize(2*p, 3*p);   
        }
    }

    send.resize(2*pid);

    indeterminate_scatterv(receive, send);

    if (send.size() == 2*pid)
    {
        ut.passes("correct number of processors in scatterv");
    }
    else
    {
        ut.failure("NOT correct number of processors in scatterv");
    }
    for (unsigned i=0; i<2*pid; ++i)
    {
        if (send[i] != 3*pid)
        {
            ut.failure("NOT correct values in scatterv");
        }
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstDeterminateGatherScatter(ut);
        tstIndeterminateGatherScatterv(ut);
        tstDeterminateGatherScatterv(ut);
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
