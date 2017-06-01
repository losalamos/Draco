//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstC4_Status.cc
 * \author Robert B. Lowrie
 * \date   Friday May 26 19:58:19 2017
 * \brief  Unit test for C4_Status class.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <vector>

using namespace std;
using rtt_c4::C4_Status;
using rtt_c4::C4_Req;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tst2Procs(rtt_dsxx::UnitTest &ut) {

  if (rtt_c4::nodes() != 2)
    return;

  PASSMSG("Running tst2Procs.");

  const int my_proc = rtt_c4::node();
  C4_Status status;
  C4_Req request;

  const int num_int = 2;
  const int num_double = 5;
  const int tag = 101;

  if (my_proc == 0) {
    std::vector<int> send_buffer(num_int);
    std::vector<double> recv_buffer(num_double);
    rtt_c4::receive_async(request, &recv_buffer[0], num_double, 1, tag);
    rtt_c4::send_async(&send_buffer[0], num_int, 1, tag);
    request.wait(&status);
    if (status.get_source() == 1)
      PASSMSG("get_source() passed on proc 0");
    else
      FAILMSG("get_source() failed on proc 0");
    if (status.get_message_size() == num_double * sizeof(double))
      PASSMSG("get_message_size() passed on proc 0");
    else
      FAILMSG("get_message_size() failed on proc 0");
    if (status.get_status_obj())
      PASSMSG("get_status_obj() passed on proc 0");
    else
      FAILMSG("get_status_obj() failed on proc 0");
  } else { // my_proc == 1
    std::vector<double> send_buffer(num_double);
    std::vector<int> recv_buffer(num_int);
    rtt_c4::receive_async(request, &recv_buffer[0], num_int, 0, tag);
    rtt_c4::send_async(&send_buffer[0], num_double, 0, tag);
    request.wait(&status);
    if (status.get_source() == 0)
      PASSMSG("get_source() passed on proc 1");
    else
      FAILMSG("get_source() failed on proc 1");
    if (status.get_message_size() == num_int * sizeof(int))
      PASSMSG("get_message_size() passed on proc 1");
    else
      FAILMSG("get_message_size() failed on proc 1");
    if (status.get_status_obj())
      PASSMSG("get_status_obj() passed on proc 1");
    else
      FAILMSG("get_status_obj() failed on proc 1");
  }

  if (ut.numFails == 0)
    PASSMSG("tstC4_Status() is okay.");

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    if (rtt_c4::nodes() == 2)
      tst2Procs(ut);
    else
      FAILMSG("tstC4_Status should only be run on 2 procs!");
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstC4_Status.cc
//---------------------------------------------------------------------------//
