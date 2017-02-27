//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstsend_is.cc
 * \author Kelly Thompson
 * \date   Friday, Dec 07, 2012, 14:02 pm
 * \brief  Unit tests for rtt_c4::send_is()
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: tstsend_is.cc 5830 2011-05-05 19:43:43Z kellyt $
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <sstream>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_simple(rtt_dsxx::UnitTest &ut) {
  // borrowed from http://mpi.deino.net/mpi_functions/MPI_Issend.html.

  if (rtt_c4::node() == 0)
    std::cout << "Test send_is() by sending data to proc myid+1..."
              << std::endl;

  // C4_Req communication handles.
  std::vector<rtt_c4::C4_Req> comm_int(2);

  // for point-to-point communiction we need to know neighbor's identifiers:
  // left, right.
  int right = (rtt_c4::node() + 1) % rtt_c4::nodes();
  int left = rtt_c4::node() - 1;
  if (left < 0)
    left = rtt_c4::nodes() - 1;

  // create some data to send/recv
  unsigned int const bsize(10);
  std::vector<int> buffer2(bsize);
  std::vector<int> buffer1(bsize);
  for (size_t i = 0; i < bsize; ++i)
    buffer1[i] = 1000 * rtt_c4::node() + i;

  // post asynchronous receives.
  comm_int[0] = rtt_c4::receive_async(&buffer2[0], buffer2.size(), left);

  try {
    // send data using non-blocking synchronous send.
    rtt_c4::send_is(comm_int[1], &buffer1[0], buffer1.size(), right);

    // wait for all communication to finish
    rtt_c4::wait_all(comm_int.size(), &comm_int[0]);

    // exected results
    std::vector<int> expected(bsize);
    for (size_t i = 0; i < bsize; ++i)
      expected[i] = 1000 * left + i;

    if (expected == buffer2) {
      std::ostringstream msg;
      msg << "Expected data found after send_is() on node " << rtt_c4::node()
          << ".";
      PASSMSG(msg.str());
    } else {
      std::ostringstream msg;
      msg << "Did not find expected data after send_is() on node "
          << rtt_c4::node() << ".";
      FAILMSG(msg.str());
    }
  } catch (rtt_dsxx::assertion const &error) {
#ifdef C4_SCALAR
    PASSMSG("Successfully caught a ds++ exception while trying to use "
            "send_is() in a C4_SCALAR build.");
#else
    FAILMSG("Encountered a ds++ exception while testing send_is().");
#endif
  }

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_simple(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstsend_is.cc
//---------------------------------------------------------------------------//
