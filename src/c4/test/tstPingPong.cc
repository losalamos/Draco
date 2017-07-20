//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstPingPong.cc
 * \author Thomas M. Evans
 * \date   Tue Apr  2 15:57:11 2002
 * \brief  Ping Pong communication test.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>

using namespace std;

using rtt_c4::C4_Req;
using rtt_c4::C4_Status;
using rtt_c4::C4_Traits;
using rtt_c4::send;
using rtt_c4::receive;
using rtt_c4::send_async;
using rtt_c4::receive_async;
using rtt_dsxx::soft_equiv;
using rtt_c4::probe;
using rtt_c4::blocking_probe;
using rtt_c4::send_receive;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void blocking_ping_pong(rtt_dsxx::UnitTest &ut) {
  if (rtt_c4::nodes() != 2)
    return;

  bool b = false;
  char c = 0;
  int i = 0;
  long l = 0;
  float f = 0;
  double d = 0;

  // assign on node 0
  if (rtt_c4::node() == 0) {
    b = true;
    c = 'A';
    i = 1;
    l = 1000;
    f = 1.5;
    d = 2.5;

    // send out data
    // Test both active and depricated forms of the send command.
    send(&b, 1, 1);
    send(&c, 1, 1);
    send(&i, 1, 1);
    send(&l, 1, 1);
    send(&f, 1, 1);
    send(&d, 1, 1);

    // receive back
    // Test both active and depricated forms of the receive command.
    receive(&b, 1, 1);
    receive(&c, 1, 1);
    receive(&i, 1, 1);
    receive(&l, 1, 1);
    receive(&f, 1, 1);
    receive(&d, 1, 1);

    // check values
    if (b != false)
      ITFAILS;
    if (c != 'B')
      ITFAILS;
    if (i != 2)
      ITFAILS;
    if (l != 2000)
      ITFAILS;
    if (!soft_equiv(f, 2.5f))
      ITFAILS;
    if (!soft_equiv(d, 3.5))
      ITFAILS;
  }

  // receive and send on node 1
  if (rtt_c4::node() == 1) {
    // receive from node 0
    receive(&b, 1, 0);
    receive(&c, 1, 0);
    receive(&i, 1, 0);
    receive(&l, 1, 0);
    receive(&f, 1, 0);
    receive(&d, 1, 0);

    // check values
    if (b != true)
      ITFAILS;
    if (c != 'A')
      ITFAILS;
    if (i != 1)
      ITFAILS;
    if (l != 1000)
      ITFAILS;
    if (!soft_equiv(f, 1.5f))
      ITFAILS;
    if (!soft_equiv(d, 2.5))
      ITFAILS;

    // assign new values
    b = false;
    c = 'B';
    i = 2;
    l = 2000;
    f = 2.5;
    d = 3.5;

    // send them back
    send(&b, 1, 0);
    send(&c, 1, 0);
    send(&i, 1, 0);
    send(&l, 1, 0);
    send(&f, 1, 0);
    send(&d, 1, 0);
  }

  rtt_c4::global_barrier();
  if (ut.numFails == 0) {
    ostringstream m;
    m << "Blocking Send/Recv communication ok on " << rtt_c4::node();
    PASSMSG(m.str());
  } else {
    ostringstream m;
    m << "Blocking Send/Recv communication failed on " << rtt_c4::node();
    FAILMSG(m.str());
  }
  return;
}
//---------------------------------------------------------------------------//

void non_blocking_ping_pong(rtt_dsxx::UnitTest &ut) {

  if (rtt_c4::nodes() != 2)
    return;
  bool b = false;
  char c = 0;
  int i = 0;
  long l = 0;
  float f = 0;
  double d = 0;

  bool br = false;
  char cr = 0;
  int ir = 0;
  long lr = 0;
  float fr = 0;
  double dr = 0;

  // send requests
  C4_Req brs, crs, irs, lrs, frs, drs;

  // receive requests
  C4_Req brr, crr, irr, lrr, frr, drr;

  // assign on node 0
  if (rtt_c4::node() == 0) {
    // post receives
    // Test two forms of the receive_async command plus one deprecated
    // form (namespace C4::)

    brr = receive_async(&br, 1, 1);
    receive_async(crr, &cr, 1, 1);
    irr = receive_async(&ir, 1, 1);
    receive_async(lrr, &lr, 1, 1);
    frr = receive_async(&fr, 1, 1);
    receive_async(drr, &dr, 1, 1);

    // give values to the send data
    b = true;
    c = 'A';
    i = 1;
    l = 1000;
    f = 1.5;
    d = 2.5;

    // send out data
    // Test two forms of the send_async command plus one deprecated
    // form (namespace C4::)
    brs = send_async(&b, 1, 1);
    send_async(crs, &c, 1, 1);
    irs = send_async(&i, 1, 1);
    send_async(lrs, &l, 1, 1);
    frs = send_async(&f, 1, 1);
    send_async(drs, &d, 1, 1);

    // wait for sends to be finished
    brs.wait();
    crs.wait();
    irs.wait();
    lrs.wait();
    frs.wait();
    drs.wait();

    // wait on receives and check

    C4_Status status;

    brr.wait(&status);
    if (status.get_message_size() != sizeof(bool))
      ITFAILS;
    if (status.get_source() != 1)
      ITFAILS;

    crr.wait(&status);
    if (status.get_message_size() != 1)
      ITFAILS;
    if (status.get_source() != 1)
      ITFAILS;

    irr.wait(&status);
    if (status.get_message_size() != sizeof(int))
      ITFAILS;
    if (status.get_source() != 1)
      ITFAILS;

    lrr.wait(&status);
    if (status.get_message_size() != sizeof(long))
      ITFAILS;
    if (status.get_source() != 1)
      ITFAILS;

    frr.wait(&status);
    if (status.get_message_size() != sizeof(float))
      ITFAILS;
    if (status.get_source() != 1)
      ITFAILS;

    drr.wait(&status);
    if (status.get_message_size() != sizeof(double))
      ITFAILS;
    if (status.get_source() != 1)
      ITFAILS;

    // check values
    if (br != false)
      ITFAILS;
    if (cr != 'B')
      ITFAILS;
    if (ir != 2)
      ITFAILS;
    if (lr != 2000)
      ITFAILS;
    if (!soft_equiv(fr, 2.5f))
      ITFAILS;
    if (!soft_equiv(dr, 3.5))
      ITFAILS;
  }

  // receive and send on node 1
  if (rtt_c4::node() == 1) {
    // post receives
    // Test both function that provide equivalent functionality.
    brr = receive_async(&br, 1, 0);
    receive_async(crr, &cr, 1, 0);
    irr = receive_async(&ir, 1, 0);
    receive_async(lrr, &lr, 1, 0);
    receive_async(frr, &fr, 1, 0);
    receive_async(drr, &dr, 1, 0);

    // check that all are inuse
    if (!brr.inuse())
      ITFAILS;
    if (!crr.inuse())
      ITFAILS;
    if (!irr.inuse())
      ITFAILS;
    if (!lrr.inuse())
      ITFAILS;
    if (!frr.inuse())
      ITFAILS;
    if (!drr.inuse())
      ITFAILS;

    // check on receives
    int done = 0;
    while (done < 6) {
      if (brr.complete())
        done++;
      if (crr.complete())
        done++;
      if (irr.complete())
        done++;
      if (lrr.complete())
        done++;
      if (frr.complete())
        done++;
      if (drr.complete())
        done++;
    }

    if (br != true)
      ITFAILS;
    if (cr != 'A')
      ITFAILS;
    if (ir != 1)
      ITFAILS;
    if (lr != 1000)
      ITFAILS;
    if (!soft_equiv(fr, 1.5f))
      ITFAILS;
    if (!soft_equiv(dr, 2.5))
      ITFAILS;

    // assign new values
    b = false;
    c = 'B';
    i = 2;
    l = 2000;
    f = 2.5;
    d = 3.5;

    // send them back
    // Test both function that provide equivalent functionality.

    brs = send_async(&b, 1, 0);
    send_async(crs, &c, 1, 0);
    irs = send_async(&i, 1, 0);
    send_async(lrs, &l, 1, 0);
    send_async(frs, &f, 1, 0);
    send_async(drs, &d, 1, 0);

    // wait for sends to be finished
    brs.wait();
    crs.wait();
    irs.wait();
    lrs.wait();
    frs.wait();
    drs.wait();
  }
  rtt_c4::global_barrier();

  // check that all requests are done
  if (brs.inuse())
    ITFAILS;
  if (crs.inuse())
    ITFAILS;
  if (irs.inuse())
    ITFAILS;
  if (lrs.inuse())
    ITFAILS;
  if (frs.inuse())
    ITFAILS;
  if (drs.inuse())
    ITFAILS;

  if (brr.inuse())
    ITFAILS;
  if (crr.inuse())
    ITFAILS;
  if (irr.inuse())
    ITFAILS;
  if (lrr.inuse())
    ITFAILS;
  if (frr.inuse())
    ITFAILS;
  if (drr.inuse())
    ITFAILS;

  if (ut.numFails == 0) {
    ostringstream m;
    m << "Non-blocking Send/Recv communication ok on " << rtt_c4::node();
    PASSMSG(m.str());
  }

  return;
}

//---------------------------------------------------------------------------//
// Test the C4_Req.free() function.
// After a asynchronous receive command is given, it can be terminated by
// using the free() command.
//---------------------------------------------------------------------------//

void tstC4_Req_free() {

  if (rtt_c4::nodes() != 2)
    return;

  // receive requests
  C4_Req crr;
  char cr = 0;

  // assign on node 0
  if (rtt_c4::node() == 0) {
    // post receives
    receive_async(crr, &cr, 1, 1);
    // void the receive request.
    crr.free();
  }

  return;
}

//---------------------------------------------------------------------------//

void probe_ping_pong(rtt_dsxx::UnitTest &ut) {

  if (rtt_c4::nodes() != 2)
    return;

  int i = 0;
  int ir = 0;

  // send requests
  C4_Req irs;

  // receive requests
  C4_Req irr;

  // Only 2 procs so, dest is either 0 or 1.
  int dest(std::abs(rtt_c4::node() - 1));

  // assign on node 1
  if (rtt_c4::node() == 1) {

    // give values to the send data
    i = 2;

    // send out data
    send_async(irs, &i, 1, dest);
  }

  // receive and send on node 0
  if (rtt_c4::node() == 0) {
    // test the probe function
    int message_size;
    for (;;) {
      if (probe(dest, C4_Traits<int *>::tag, message_size)) {
        if (message_size == sizeof(int)) {
          PASSMSG("Probe returned correct size");
        } else {
          FAILMSG("Probe returned WRONG size");
        }
        break;
      }
    }

    // test the blocking probe function
    blocking_probe(dest, C4_Traits<int *>::tag, message_size);
    if (message_size == sizeof(int)) {
      PASSMSG("Blocking probe returned correct size");
    } else {
      FAILMSG("Blocking probe returned WRONG size");
    }

    // post receives
    receive_async(irr, &ir, 1, dest);
  }

  rtt_c4::global_barrier();

  if (ut.numFails == 0) {
    ostringstream m;
    m << "Probe communication ok on " << rtt_c4::node();
    PASSMSG(m.str());
  }
  return;
}

void send_receive_ping_pong(rtt_dsxx::UnitTest &ut) {
  if (rtt_c4::nodes() != 2)
    return;

  bool b = false, br;
  char c = 0, cr;
  int i = 0, ir;
  long l = 0, lr;
  float f = 0, fr;
  double d = 0, dr;

  // assign on node 0
  if (rtt_c4::node() == 0) {
    b = true;
    c = 'A';
    i = 1;
    l = 1000;
    f = 1.5;
    d = 2.5;

    send_receive(&b, 1, 1, &br, 1, 1);
    send_receive(&c, 1, 1, &cr, 1, 1);
    send_receive(&i, 1, 1, &ir, 1, 1);
    send_receive(&l, 1, 1, &lr, 1, 1);
    send_receive(&f, 1, 1, &fr, 1, 1);
    send_receive(&d, 1, 1, &dr, 1, 1);

    // check values
    if (br != false)
      ITFAILS;
    if (cr != 'B')
      ITFAILS;
    if (ir != 2)
      ITFAILS;
    if (lr != 2000)
      ITFAILS;
    if (!soft_equiv(fr, 2.5f))
      ITFAILS;
    if (!soft_equiv(dr, 3.5))
      ITFAILS;
  }

  // receive and send on node 1
  if (rtt_c4::node() == 1) {
    // assign new values
    b = false;
    c = 'B';
    i = 2;
    l = 2000;
    f = 2.5;
    d = 3.5;

    send_receive(&b, 1, 0, &br, 1, 0);
    send_receive(&c, 1, 0, &cr, 1, 0);
    send_receive(&i, 1, 0, &ir, 1, 0);
    send_receive(&l, 1, 0, &lr, 1, 0);
    send_receive(&f, 1, 0, &fr, 1, 0);
    send_receive(&d, 1, 0, &dr, 1, 0);

    // check values
    if (br != true)
      ITFAILS;
    if (cr != 'A')
      ITFAILS;
    if (ir != 1)
      ITFAILS;
    if (lr != 1000)
      ITFAILS;
    if (!soft_equiv(fr, 1.5f))
      ITFAILS;
    if (!soft_equiv(dr, 2.5))
      ITFAILS;
  }

  rtt_c4::global_barrier();
  if (ut.numFails == 0) {
    ostringstream m;
    m << "send-receive communication ok on " << rtt_c4::node();
    PASSMSG(m.str());
  } else {
    ostringstream m;
    m << "send-receive communication failed on " << rtt_c4::node();
    FAILMSG(m.str());
  }
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    std::cout << "This is " << rtt_c4::get_processor_name() << std::endl;

    Insist(rtt_c4::nodes() == 2,
           "This designed is designed for exactly 2 Processors.");

    blocking_ping_pong(ut);
    non_blocking_ping_pong(ut);
    probe_ping_pong(ut);
    send_receive_ping_pong(ut);
    tstC4_Req_free();
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstPingPong.cc
//---------------------------------------------------------------------------//
