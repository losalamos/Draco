//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstTigsTrace.hh
 * \author Bob Webster, Thomas M. Evans
 * \date   Fri Sep 30 12:53:25 2005
 * \brief  Test the TigsTrace class.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * These comm classes are used to perform gather/scatter operations.  They
 * were originally developed by Bob Webster for the Antero Project.  They have
 * been imported here and (slightly) refactored for use in wedgehog
 * census-mapping operations.  They are used with permission of the original
 * author.
 */
//---------------------------------------------------------------------------//
// $Id: Tigs.hh 6056 2012-06-19 19:05:27Z kellyt $
//---------------------------------------------------------------------------//

#include "c4/TigsTrace.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <functional>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void test_trace( rtt_dsxx::UnitTest &ut )
{
    // bool debug=false;

    size_t numProcs=rtt_c4::nodes();

    // ------------------------------------------------------------------------
    // First step is to initialize a trace for a fairly complex communication
    // pattern. The processor to processor communication is as follows...
    //
    //   n-1 -> 0
    //   n-2 -> 1, 0
    //   n-3 -> 2, 1, 0
    //   n-4 -> 3, 2, 1, 0
    //       .
    //       .
    //       .
    //   0   -> n-1, ... , 3, 2, 1, 0
    //
    //  Additionally, there will be procRank()+numProcs() elements of the
    //  range data on each processor.
    //
    size_t const nRange=rtt_c4::node()+rtt_c4::nodes();
    std::vector<int> JstarTest(rtt_c4::nodes()+1,0);
    size_t runningTot(0);
    for( int pe=0; pe<rtt_c4::nodes()+1; pe++ )
    {
        JstarTest[pe]=runningTot;
        runningTot+=rtt_c4::nodes()+pe;
    }

    // reverse the pattern for the communication domain data
    size_t inversePE=rtt_c4::nodes()-rtt_c4::node();
    size_t I=JstarTest[inversePE];
    std::vector<int> M(I, 0);
    for( size_t i=0; i<I; i++ )
        M[i]=i;

    rtt_c4::global_barrier();
    rtt_c4::TigsTrace itrace(M, nRange);

    std::vector<int> Jlocal(nRange,0);
    int base=0;
    int offsetTest=10;
    for( int pe=0; pe<rtt_c4::node(); pe++)
        base+=pe+rtt_c4::nodes();
    for( size_t i=0; i<nRange; i++)
        Jlocal[i]=base+i + offsetTest;

    //----------------------------------------------------------------------
    // test gather

    rtt_c4::global_barrier();
    if( rtt_c4::node()==0 )
        std::cout << "Gather Test" << std::endl;

    std::vector<int> JGather(M.size(),0);

    itrace.gather( Jlocal.begin(), Jlocal.end(),
                   JGather.begin(), JGather.end() );

    int ierrGather=0;
    for( size_t i=0; i<M.size(); i++ )
        if( M[i] != (JGather[i]-offsetTest) )
            ierrGather=1;

    if( ierrGather==0 )
        PASSMSG("gather std::vector<int> works.");
    else
        FAILMSG("gather std::vector<int> failed to return expected array.");

    //----------------------------------------------------------------------
    //  test the scatterList function

    rtt_c4::global_barrier();
    if( rtt_c4::node()==0 )
        std::cout << "ScatterList Test" << std::endl;

    int ierrScatterList=0;
    size_t rangeSize=nRange;
    size_t listSize=itrace.getListSize();
    // int rangeSize=itrace.getRangeSize();
    std::vector<int> Ilocal(I,0);
    std::vector<int> IScatter(listSize,0);
    std::vector<int> counts_ret(rangeSize);;

    // Enumerate the i locations so that when we scatter we know where they
    // went.
    int baseI=0;
    for( int pe=0; pe<rtt_c4::node() ; pe++ )
    {
        baseI += JstarTest[rtt_c4::nodes()-pe];
    }
    for( size_t i=0; i<I; i++ )
        Ilocal[i]=1000*rtt_c4::node()+i+1000;

    itrace.scatterList(Ilocal.begin(),     Ilocal.end(),
                       counts_ret.begin(), counts_ret.end(),
                       IScatter.begin(),   IScatter.end());

    // Test the result:
    for( size_t pe=0; pe<numProcs; pe++ )
    {
        if( rtt_c4::node() == static_cast<int>(pe) )
        {
            int k=0;
            for( size_t j=0; j<counts_ret.size(); j++ )
            {
                for( int l=0; l<counts_ret[j]; l++ )
                {
                    if( IScatter[k+l] != static_cast<int>(1000*(l+1)
                                                          + j +base) )
                        ierrScatterList=1;
                }
                k+=counts_ret[j];
            }
        }
        rtt_c4::global_barrier();
    }

    if( ierrScatterList==0 )
        PASSMSG( "scatterList std::vector<int> works." );
    else
        FAILMSG( "scatterList std::vector<int> failed to return expected array." );

    //----------------------------------------------------------------------
    //  test a scatter sum

    rtt_c4::global_barrier();
    if( rtt_c4::node()==0 )
        std::cout << "Scatter Test" << std::endl;

    int ierrScatter=0;
    std::vector<int> ICount(I, 1);
    std::vector<int> SSumResult(nRange,0);

    itrace.scatter( ICount.begin(),   ICount.end(),
                    SSumResult.begin(), SSumResult.end(),
                    std::plus<int>() );

    // Test the result:
    for( size_t pe=0; pe<numProcs; pe++ )
    {
        if( rtt_c4::node() == static_cast<int>(pe) )
            for( size_t j=0; j<counts_ret.size(); j++ )
                if( (SSumResult[j]+rtt_c4::node()) != rtt_c4::nodes() )
                    ierrScatter=1;
        rtt_c4::global_barrier();
    }

    if( ierrScatter==0 )
        PASSMSG( "scatter std::vector<int> works." );
    else
        FAILMSG( "scatter std::vector<int> failed to return expected array." );

    return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
    try { test_trace(ut); }
    UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstTigsTrace.cc
//---------------------------------------------------------------------------//
