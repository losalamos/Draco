//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstBroadcast.cc
 * \author Thomas M. Evans
 * \date   Tue Apr  2 15:57:11 2002
 * \brief  Ping Pong communication test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4_test.hh"
#include "../C4_Traits.hh"
#include "../Release.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <sstream>

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
    char   c = 0;
    int    i = 0;
    long   l = 0;
    float  f = 0;
    double d = 0;

    // assign on node 0
    if (rtt_c4::node() == 0)
    {
	c = 'A';
	i = 1;
	l = 1000;
	f = 1.5;
	d = 2.5;
    }

    // send out data, using node 0 as root
    broadcast(&c, 1, 0);
    broadcast(&i, 1, 0);
    broadcast(&l, 1, 0);
    broadcast(&f, 1, 0);
    broadcast(&d, 1, 0);

    // check values
    if (c != 'A')             ITFAILS;
    if (i != 1)               ITFAILS;
    if (l != 1000)            ITFAILS;
    if (!soft_equiv(f, 1.5f)) ITFAILS;
    if (!soft_equiv(d, 2.5))  ITFAILS;

    rtt_c4::global_barrier();

    if (rtt_c4_test::passed)
    {
	ostringstream m;
	m << "Broadcast communication ok on " << rtt_c4::node();
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
		cout << argv[0] << ": version " << rtt_c4::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	test_simple();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstBroadcast, " << ass.what()
	     << endl;
	rtt_c4::finalize();
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
