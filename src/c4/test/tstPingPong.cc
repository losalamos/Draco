//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstPingPong.cc
 * \author Thomas M. Evans
 * \date   Tue Apr  2 15:57:11 2002
 * \brief  Ping Pong communication test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4_test.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

using rtt_c4::C4_Req;
using rtt_c4::send;
using rtt_c4::receive;
using rtt_c4::send_async;
using rtt_c4::receive_async;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void blocking_ping_pong()
{
    if (rtt_c4::nodes() != 2) return;
    
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
	
	// send out data
	send(&c, 1, 1);
	send(&i, 1, 1);
	send(&l, 1, 1);
	send(&f, 1, 1);
	send(&d, 1, 1);

	// receive back
	receive(&c, 1, 1);
	receive(&i, 1, 1);
	receive(&l, 1, 1);
	receive(&f, 1, 1);
	receive(&d, 1, 1);

	// check values
	if (c != 'B')             ITFAILS;
	if (i != 2)               ITFAILS;
	if (l != 2000)            ITFAILS;
	if (!soft_equiv(f, 2.5f)) ITFAILS;
	if (!soft_equiv(d, 3.5))  ITFAILS;
    }

    // receive and send on node 1
    if (rtt_c4::node() == 1)
    {
	// receive from node 0
	receive(&c, 1, 0);
	receive(&i, 1, 0);
	receive(&l, 1, 0);
	receive(&f, 1, 0);
	receive(&d, 1, 0);

	// check values
	if (c != 'A')             ITFAILS;
	if (i != 1)               ITFAILS;
	if (l != 1000)            ITFAILS;
	if (!soft_equiv(f, 1.5f)) ITFAILS;
	if (!soft_equiv(d, 2.5))  ITFAILS;

	// assign new values
	c = 'B';
	i = 2;
	l = 2000;
	f = 2.5;
	d = 3.5;

	// send them back
	send(&c, 1, 0);
	send(&i, 1, 0);
	send(&l, 1, 0);
	send(&f, 1, 0);
	send(&d, 1, 0);
    }

    rtt_c4::global_barrier();

    if (rtt_c4_test::passed)
    {
	ostringstream m;
	m << "Blocking Send/Recv communication ok on " << rtt_c4::node();
	PASSMSG(m.str());
    }
}
//---------------------------------------------------------------------------//

void non_blocking_ping_pong()
{
    if (rtt_c4::nodes() != 2) return;
    
    char   c = 0;
    int    i = 0;
    long   l = 0;
    float  f = 0;
    double d = 0;
    
    char   cr = 0;
    int    ir = 0;
    long   lr = 0;
    float  fr = 0;
    double dr = 0;

    // send requests
    C4_Req crs, irs, lrs, frs, drs;

    // receive requests
    C4_Req crr, irr, lrr, frr, drr;

    // assign on node 0
    if (rtt_c4::node() == 0)
    {
	// post receives
	receive_async(crr, &cr, 1, 1);
	receive_async(irr, &ir, 1, 1);
	receive_async(lrr, &lr, 1, 1);
	receive_async(frr, &fr, 1, 1);
	receive_async(drr, &dr, 1, 1);

	// give values to the send data
	c = 'A';
	i = 1;
	l = 1000;
	f = 1.5;
	d = 2.5;
	
	// send out data
	send_async(crs, &c, 1, 1);
	send_async(irs, &i, 1, 1);
	send_async(lrs, &l, 1, 1);
	send_async(frs, &f, 1, 1);
	send_async(drs, &d, 1, 1);

	// wait for sends to be finished
	crs.wait();
	irs.wait();
	lrs.wait();
	frs.wait();
	drs.wait();

	// wait on receives and check
	crr.wait();
	irr.wait();
	lrr.wait();
	frr.wait();
	drr.wait();

	// check values
	if (cr != 'B')             ITFAILS;
	if (ir != 2)               ITFAILS;
	if (lr != 2000)            ITFAILS;
	if (!soft_equiv(fr, 2.5f)) ITFAILS;
	if (!soft_equiv(dr, 3.5))  ITFAILS;
    }

    // receive and send on node 1
    if (rtt_c4::node() == 1)
    {
	// post receives
	receive_async(crr, &cr, 1, 0);
	receive_async(irr, &ir, 1, 0);
	receive_async(lrr, &lr, 1, 0);
	receive_async(frr, &fr, 1, 0);
	receive_async(drr, &dr, 1, 0);

	// check that all are inuse
	if (!crr.inuse()) ITFAILS;
	if (!irr.inuse()) ITFAILS;
	if (!lrr.inuse()) ITFAILS;
	if (!frr.inuse()) ITFAILS;
	if (!drr.inuse()) ITFAILS;

	// check on receives
	int done = 0;
	while (done < 5)
	{
	    if (crr.complete()) done++;
	    if (irr.complete()) done++;
	    if (lrr.complete()) done++;
	    if (frr.complete()) done++;
	    if (drr.complete()) done++;
	}

	if (cr != 'A')             ITFAILS;
	if (ir != 1)               ITFAILS;
	if (lr != 1000)            ITFAILS;
	if (!soft_equiv(fr, 1.5f)) ITFAILS;
	if (!soft_equiv(dr, 2.5))  ITFAILS;

	// assign new values
	c = 'B';
	i = 2;
	l = 2000;
	f = 2.5;
	d = 3.5;

	// send them back
	send_async(crs, &c, 1, 0);
	send_async(irs, &i, 1, 0);
	send_async(lrs, &l, 1, 0);
	send_async(frs, &f, 1, 0);
	send_async(drs, &d, 1, 0);

	// wait for sends to be finished
	crs.wait();
	irs.wait();
	lrs.wait();
	frs.wait();
	drs.wait();
    }

    rtt_c4::global_barrier();

    // check that all requests are done
    if (crs.inuse()) ITFAILS;
    if (irs.inuse()) ITFAILS;
    if (lrs.inuse()) ITFAILS;
    if (frs.inuse()) ITFAILS;
    if (drs.inuse()) ITFAILS;

    if (crr.inuse()) ITFAILS;
    if (irr.inuse()) ITFAILS;
    if (lrr.inuse()) ITFAILS;
    if (frr.inuse()) ITFAILS;
    if (drr.inuse()) ITFAILS;

    if (rtt_c4_test::passed)
    {
	ostringstream m;
	m << "Non-blocking Send/Recv communication ok on " << rtt_c4::node();
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
	blocking_ping_pong();

	non_blocking_ping_pong();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstPingPong, " << ass.what()
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
	    cout << "**** tstPingPong Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstPingPong on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstPingPong.cc
//---------------------------------------------------------------------------//
