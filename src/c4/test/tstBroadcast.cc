//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstBroadcast.cc
 * \author Thomas M. Evans
 * \date   Tue Apr  2 15:57:11 2002
 * \brief  Ping Pong communication test.
 * \note   Copyright (C) 2002-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../C4_Traits.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "../ParallelUnitTest.hh"

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"

#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

using rtt_c4::C4_Req;
using rtt_c4::C4_Traits;
using rtt_c4::broadcast;
using rtt_dsxx::soft_equiv;

#define PASSMSG(A) ut.passes(A)
#define FAILMSG(A) ut.failure(A)
#define ITFAILS    ut.failure( __LINE__ )

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_simple( rtt_dsxx::UnitTest &ut )
{
    using std::vector;
    
    char   c = 0;
    int    i = 0;
    long   l = 0;
    float  f = 0;
    double d = 0;
    vector<double> vref(10,3.1415);
    vector<double> v(10,0.0);
    
    // assign on node 0
    if (rtt_c4::node() == 0)
    {
	c = 'A';
	i = 1;
	l = 1000;
	f = 1.5;
	d = 2.5;
        v = vref;
    }
    
    // send out data, using node 0 as root
    broadcast(&c, 1, 0);
    broadcast(&i, 1, 0);
    broadcast(&l, 1, 0);
    broadcast(&f, 1, 0);
    broadcast(&d, 1, 0);
    broadcast(v.begin(),v.end(),v.begin());

    // check values
    if (c != 'A')             ITFAILS;
    if (i != 1)               ITFAILS;
    if (l != 1000)            ITFAILS;
    if (!soft_equiv(f, 1.5f)) ITFAILS;
    if (!soft_equiv(d, 2.5))  ITFAILS;
    if (!soft_equiv(v.begin(),v.end(),vref.begin(),vref.end()))  ITFAILS;
    
    rtt_c4::global_barrier();

    if( ut.numFails == 0 )
    {
	ostringstream m;
	m << "test_simple() ok on " << rtt_c4::node();
	PASSMSG(m.str());
    }
    return;
}

//---------------------------------------------------------------------------//
//  By adjusting the parameters below, this test will overflow the MPI memory
//  buffers.  Read the comments below if you'd like to do this.
void test_loop( rtt_dsxx::UnitTest &ut )
{
    // save state
    unsigned const nf (ut.numFails);
    
    // >>> kmax controls how much data is broadcast.  If kmax is too big
    // >>> (like 10000000), shmem will fail.
    int const kmax = 10;
    
    if (rtt_c4::node() == 0) // host proc
    {
	// send out the values on host
	for ( int k = 0; k < kmax; ++k )
	{
	    Insist(! broadcast(&k, 1, 0), "MPI Error");
	    double foo = k + 0.5;
	    Insist(! broadcast(&foo, 1, 0), "MPI Error");
	}
    }
    else // all other procs
    {
	// >>> Use sleep() if you want the host processor to fill up the
	// >>> buffers.  We comment out the sleep() command here because it's
	// >>> not supported on all all platforms.

	// sleep(10);
	
	int kk;
	double foofoo;
	for ( int k = 0; k < kmax; ++k )
	{
	    kk = -1;
	    foofoo = -2.0;
	    Insist(! broadcast(&kk, 1, 0), "MPI Error");
	    if ( kk != k ) ITFAILS;
	    Insist(! broadcast(&foofoo, 1, 0), "MPI Error");
	    if ( foofoo != k + 0.5 ) ITFAILS;
	}
    }

    if( ut.numFails == nf )
    {
	ostringstream m;
	m << "test_loop() ok on " << rtt_c4::node();
	PASSMSG(m.str());
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        // Unit tests
        test_simple(ut);
        test_loop(  ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstBroadcast, " 
                  << err.what() << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstBroadcast, " 
                  << "An unknown exception was thrown." << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstBroadcast.cc
//---------------------------------------------------------------------------//
