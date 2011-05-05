//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstBroadcast.cc
 * \author Thomas M. Evans
 * \date   Tue Apr  2 15:57:11 2002
 * \brief  Ping Pong communication test.
 * \note   Copyright (C) 2006-2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <sstream>

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include "c4_test.hh"
#include "../C4_Traits.hh"
#include "ds++/Release.hh"
#include "../global.hh"
#include "../SpinLock.hh"

using namespace std;

using rtt_c4::C4_Req;
using rtt_c4::C4_Traits;
using rtt_c4::broadcast;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_simple()
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

    if (rtt_c4_test::passed)
    {
	ostringstream m;
	m << "test_simple() ok on " << rtt_c4::node();
	PASSMSG(m.str());
    }
}

//---------------------------------------------------------------------------//
//  By adjusting the parameters below, this test will overflow the MPI memory
//  buffers.  Read the comments below if you'd like to do this.
void test_loop()
{
    // >>> kmax controls how much data is broadcast.  If kmax is too big
    // >>> (like 10000000), shmem will fail.
    const int kmax = 10;
    
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

    if (rtt_c4_test::passed)
    {
	ostringstream m;
	m << "test_loop() ok on " << rtt_c4::node();
	PASSMSG(m.str());
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " << rtt_dsxx::release()
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	test_simple();
	test_loop();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstBroadcast, " << ass.what()
	     << endl;
	rtt_c4::abort();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_c4_test::passed) 
	{
	    cout << "**** tstBroadcast Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstBroadcast on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstBroadcast.cc
//---------------------------------------------------------------------------//
