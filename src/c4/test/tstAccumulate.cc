//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstAccumulate.cc
 * \author Kelly Thompson
 * \date   Thursday, Nov 01, 2012, 09:14 am
 * \brief  Test c4::accumulatev function
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: tstAccumulate.cc 6288 2011-12-04 03:43:52Z kellyt $
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/accumulatev.hh"
#include "c4/Functors.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstAccumulate(UnitTest &ut)
{
    unsigned pid = node();
    unsigned const number_of_processors = nodes();

    {
        // create an array that is unique to each PE.
        vector<int> pids(number_of_processors,-1);
        pids[pid] = 1000+pid;

        // Accumulate 'max' data across PEs
        int initval(10);
        accumulatev( pids.begin(), pids.end(), initval, rtt_c4::max<int>() );
        
        std::vector<int> expected(number_of_processors);
        for( size_t i=0; i<number_of_processors; ++i )
            expected[i] = 1000+i;
        
        if( expected == pids )
            PASSMSG("accumulatev<int,max> works.");
        else
            FAILMSG("acumulatev<int,max> fails.");
    }

    {
        // create an array that is unique to each PE.
        vector<double> data(number_of_processors,-1);
        data[pid] = 1000.0+pid;

        // Accumulate 'min' data across PEs
        double initval(10000.0);
        accumulatev( data.begin(), data.end(), initval, rtt_c4::min<double>() );
        
        std::vector<double> expected(number_of_processors,
                                     number_of_processors>1?-1.0:1000.0);
        
        if( rtt_dsxx::soft_equiv(expected.begin(),expected.end(),
                                 data.begin(), data.end()) )
            PASSMSG("accumulatev<double,min> works.");
        else
            FAILMSG("acumulatev<dobule,min> fails.");
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try { tstAccumulate(ut); }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstAccumulate.cc
//---------------------------------------------------------------------------//
