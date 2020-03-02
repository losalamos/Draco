//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/test/tstsend_is.cc
 * \author Kelly Thompson
 * \date   Friday, Dec 07, 2012, 14:02 pm
 * \brief  Unit tests for rtt_c4::send_is()
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <numeric>
#include <sstream>

using rtt_dsxx::soft_equiv;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

// This is a simple class that has a static MPI type and a method to commit that
// type
class Custom {

public:
  Custom(int rank) {
    my_ints[0] = rank;
    my_ints[1] = rank * 1000;
    my_ints[2] = rank * 10000;
    my_doubles[0] = static_cast<double>(rank);
    my_doubles[1] = static_cast<double>(rank * 1000);
    my_longs[0] = rank + 100000000l;
    my_longs[1] = rank + 1000000000l;
  }
  ~Custom() {}

public:
  static const int mpi_tag = 512;
#ifdef C4_SCALAR
  static void commit_mpi_type(void) {}
  static int MPI_Type;
#else
  static MPI_Datatype MPI_Type;
  static void commit_mpi_type(void) {
    MPI_Datatype og_MPI_Custom;

    const int custom_entry_count = 3;

    // set the number of entries for each datatype
    int num_int(4);
    int num_double(2);
    int num_long(2);
    int custom_array_of_block_length[3] = {num_int, num_double, num_long};

    int int_size, double_size;
    MPI_Type_size(MPI_INT, &int_size);
    MPI_Type_size(MPI_DOUBLE, &double_size);

    // Displacements of each type in the cell
    MPI_Aint custom_array_of_block_displace[3] = {
        0, num_int * int_size, num_int * int_size + num_double * double_size};

    //Type of each memory block
    MPI_Datatype custom_array_of_types[3] = {MPI_INT, MPI_DOUBLE, MPI_LONG};

    MPI_Type_create_struct(custom_entry_count, custom_array_of_block_length,
                           custom_array_of_block_displace,
                           custom_array_of_types, &og_MPI_Custom);

    // Commit the type to MPI so it recognizes it in communication calls
    MPI_Type_commit(&og_MPI_Custom);

    // Duplicate the type so it's recognized when returned out of this
    // context (I don't know why this is necessary)
    MPI_Type_dup(og_MPI_Custom, &MPI_Type);
  }
#endif

  int get_int1(void) const { return my_ints[0]; }
  int get_int2(void) const { return my_ints[1]; }
  int get_int3(void) const { return my_ints[2]; }
  double get_double1(void) const { return my_doubles[0]; }
  double get_double2(void) const { return my_doubles[1]; }
  long get_long1(void) const { return my_longs[0]; }
  long get_long2(void) const { return my_longs[1]; }

private:
  int my_ints[3];
  double my_doubles[2];
  long my_longs[2];
};

#ifdef C4_SCALAR
int Custom::MPI_Type = 0;
#else
// the static data member needs to be defined outside the class
MPI_Datatype Custom::MPI_Type = MPI_Datatype();
#endif

//----------------------------------------------------------------------------//
void test_simple(rtt_dsxx::UnitTest &ut) {
  // borrowed from http://mpi.deino.net/mpi_functions/MPI_Issend.html.
  using namespace std;
  int const pid = rtt_c4::node();

  if (pid == 0)
    cout << "Test send_is() by sending data to proc myid+1..." << endl;

  // for point-to-point communiction we need to know neighbor's identifiers:
  // left, right.
  int right = (pid + 1) % rtt_c4::nodes();
  int left = pid - 1;
  if (left < 0)
    left = rtt_c4::nodes() - 1;
  unsigned const bsize(10);

  { // T = int

    if (pid == 0)
      std::cout << "\nStarting send_is<int> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<int> buffer2(bsize);
    vector<int> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = 1000 * pid + static_cast<int>(i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);

    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<int> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = 1000 * left + static_cast<int>(i);
      }

      if (expected == buffer2) {
        ostringstream msg;
        msg << "Expected int data found after send_is() on node " << pid << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected int data after send_is() on node " << pid
            << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<int>() in a C4_SCALAR build.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<int>().");
#endif
    }
  }

  { // T=double

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<double> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<double> buffer2(bsize);
    vector<double> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<double>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<double> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<double>(1000 * left + i);
      }

      if (soft_equiv(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected double data found after send_is() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected double data after send_is() on node "
            << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<double>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<double>().");
#endif
    }
  }

  { // T=float

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<float> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<float> buffer2(bsize);
    vector<float> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<float>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<float> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<float>(1000 * left + i);
      }

      float const eps = 1.0e-6f;
      if (soft_equiv(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end(), eps)) {
        ostringstream msg;
        msg << "Expected float data found after send_is() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected float data after send_is() on node "
            << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<float>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<float>().");
#endif
    }
  }

  { // T=long double

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<long double> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<long double> buffer2(bsize);
    vector<long double> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<long double>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<long double> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<long double>(1000 * left + i);
      }

      long double const eps = 1.0e-6f;
      if (soft_equiv(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end(), eps)) {
        ostringstream msg;
        msg << "Expected long double data found after send_is() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected long double data after send_is() on node "
            << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<long double>() in a C4_SCALAR buildb.");
#else
      FAILMSG(
          "Encountered a ds++ exception while testing send_is<long double>().");
#endif
    }
  }
  { // T=unsigned int

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<unsigned int> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<unsigned int> buffer2(bsize);
    vector<unsigned int> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<unsigned int>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<unsigned int> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<unsigned int>(1000 * left + i);
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected unsigned int data found after send_is() on node "
            << pid << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected unsigned int data after send_is() on "
               "node "
            << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<unsigned int>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "int>().");
#endif
    }
  }

  { // T=unsigned long

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<unsigned long> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<unsigned long> buffer2(bsize);
    vector<unsigned long> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<unsigned long>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<unsigned long> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<unsigned long>(1000 * left + i);
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected unsigned long data found after send_is() on node "
            << pid << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected unsigned long data after send_is() on "
               "node "
            << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<unsigned long>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "long>().");
#endif
    }
  }

  { // T=unsigned short

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<unsigned short> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<unsigned short> buffer2(bsize);
    vector<unsigned short> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<unsigned short>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<unsigned short> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<unsigned short>(1000 * left + i);
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected unsigned short data found after send_is() on node "
            << pid << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected unsigned short data after send_is() on "
               "node "
            << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<unsigned short>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "long>().");
#endif
    }
  }

  { // T=unsigned long long

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<unsigned long long> tests..."
                << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<unsigned long long> buffer2(bsize);
    vector<unsigned long long> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<unsigned long long>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<unsigned long long> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<unsigned long long>(1000 * left + i);
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected unsigned long long data found after send_is() on node "
            << pid << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected unsigned long long data after send_is() "
            << "on node " << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<unsigned long long>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "long>().");
#endif
    }
  }
  { // T=long

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<long> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<long> buffer2(bsize);
    vector<long> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<long>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<long> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<long>(1000 * left + i);
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected long data found after send_is() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected long data after send_is() "
            << "on node " << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<long>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "long>().");
#endif
    }
  }
  { // T=short

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<short> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<short> buffer2(bsize);
    vector<short> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<short>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<short> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<short>(1000 * left + i);
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected short data found after send_is() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected short data after send_is() "
            << "on node " << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<short>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "short>().");
#endif
    }
  }

  { // T=long long

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<long long> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<long long> buffer2(bsize);
    vector<long long> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = static_cast<long long>(1000 * pid + i);
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<long long> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = static_cast<long long>(1000 * left + i);
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected long long data found after send_is() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected long long data after send_is() "
            << "on node " << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<long long>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "long long>().");
#endif
    }
  }

  { // T=bool

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<bool> tests..." << std::endl;

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    bool buffer2[10];
    bool buffer1[10];
    for (unsigned i = 0; i < 10; ++i) {
      buffer1[i] = i > 5;
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0], 10, left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], 10, right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<bool> expected(10);
      for (size_t i = 0; i < 10; ++i) {
        expected[i] = i > 5;
      }

      if (std::equal(begin(expected), end(expected), begin(buffer2),
                     end(buffer2))) {
        ostringstream msg;
        msg << "Expected bool data found after send_is<bool>() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected bool data after send_is<bool>() "
            << "on node " << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<bool>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<unsigned "
              "bool>().");
#endif
    }
  }

  { // T=char

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<char> tests..." << std::endl;

    vector<char> alphabet(bsize + rtt_c4::nranks());
    std::iota(alphabet.begin(), alphabet.end(), 'A');

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<char> buffer2(bsize);
    vector<char> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = alphabet[pid + i];
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<char> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = alphabet[left + i];
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected char data found after send_is<char>() on node " << pid
            << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected char data after send_is<char>() "
            << "on node " << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<char>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<"
              "char>().");
#endif
    }
  }

  { // T=unsigned char

    rtt_c4::global_barrier();
    if (pid == 0)
      std::cout << "\nStarting send_is<unsigned char> tests..." << std::endl;

    vector<unsigned char> alphabet(bsize + rtt_c4::nranks());
    std::iota(alphabet.begin(), alphabet.end(), 'A');

    // C4_Req communication handles.
    vector<rtt_c4::C4_Req> comm(2);

    // create some size 10 data to send/recv.
    vector<unsigned char> buffer2(bsize);
    vector<unsigned char> buffer1(bsize);
    for (unsigned i = 0; i < bsize; ++i) {
      buffer1[i] = alphabet[pid + i];
    }

    // post asynchronous receives.
    comm[0] = rtt_c4::receive_async(&buffer2[0],
                                    static_cast<int>(buffer2.size()), left);
    try {
      // send data using non-blocking synchronous send.
      rtt_c4::send_is(comm[1], &buffer1[0], static_cast<int>(buffer1.size()),
                      right);

      // wait for all communication to finish
      rtt_c4::wait_all(static_cast<unsigned>(comm.size()), &comm[0]);

      // exected results
      vector<unsigned char> expected(bsize);
      for (size_t i = 0; i < bsize; ++i) {
        expected[i] = alphabet[left + i];
      }

      if (std::equal(expected.begin(), expected.end(), buffer2.begin(),
                     buffer2.end())) {
        ostringstream msg;
        msg << "Expected unsigned char data found after send_is<unsigned "
            << "char>() on node " << pid << ".";
        PASSMSG(msg.str());
      } else {
        ostringstream msg;
        msg << "Did not find expected unsigned char data after "
               "send_is<unsigned char>() on node "
            << pid << ".";
        FAILMSG(msg.str());
      }
    } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
      PASSMSG("Successfully caught a ds++ exception while trying to use "
              "send_is<unsigned char>() in a C4_SCALAR buildb.");
#else
      FAILMSG("Encountered a ds++ exception while testing send_is<"
              "unsigned char>().");
#endif
    }
  }

  rtt_c4::global_barrier();
  return;
}

//----------------------------------------------------------------------------//
void test_send_custom(rtt_dsxx::UnitTest &ut) {
  // borrowed from http://mpi.deino.net/mpi_functions/MPI_Issend.html.

  // commit the MPI type for the Custom class. This must be done before
  // send_is_custom is called. DMC checks will throw if the type has not been
  // committed because size comparison will fail and MPI throws an error when an
  // uncommited type is used in a send/receive
  Custom::commit_mpi_type();

#ifdef C4_SCALAR
  std::cout << "\nSERIAL TEST: This should just end in a try block catching ";
  std::cout << "an insist failure" << std::endl;
#else
  if (rtt_c4::node() == 0) {
    std::cout << "\nTest send_is_custom() by sending data to proc myid+1..."
              << std::endl;
    int custom_mpi_type_size(0);
    MPI_Type_size(Custom::MPI_Type, &custom_mpi_type_size);
    std::cout << " Size of custom type: " << sizeof(Custom) << std::endl;
    std::cout << " Size of custom MPI type: " << custom_mpi_type_size
              << std::endl;

    FAIL_IF_NOT(custom_mpi_type_size == sizeof(Custom));
  }
#endif

  // C4_Req communication handles.
  std::vector<rtt_c4::C4_Req> comm_int(2);

  // for point-to-point communiction we need to know neighbor's identifiers:
  // left, right.
  int right = (rtt_c4::node() + 1) % rtt_c4::nodes();
  int left = rtt_c4::node() - 1;
  if (left < 0)
    left = rtt_c4::nodes() - 1;

  // create some data to send/recv
  Custom my_custom_object(rtt_c4::node());

  // post asynchronous receives.
  Custom recv_custom_object(-1);
  rtt_c4::receive_async_custom(comm_int[0], &recv_custom_object, 1, left,
                               Custom::mpi_tag);

  try {
    // send data using non-blocking synchronous send. Custom sends check to make
    // sure that the type, T is the same size as its MPI type
    rtt_c4::send_is_custom(comm_int[1], &my_custom_object, 1, right,
                           Custom::mpi_tag);

    // make status object to get the size of the received buffer
    rtt_c4::C4_Status recv_custom_status;

    // wait for all communication to finish
    comm_int[1].wait();
    comm_int[0].wait(&recv_custom_status);

    // get the size of the message (number of objects sent) using a C4_Status
    int recv_size =
        rtt_c4::message_size_custom(recv_custom_status, Custom::MPI_Type);

    // make sure only one object was received
    FAIL_IF_NOT(recv_size == 1);

    // check that the exected results match the custom type from the left rank
    Custom expected_custom(left);

    std::cout << "Expected ints: " << expected_custom.get_int1() << " "
              << expected_custom.get_int2() << " " << expected_custom.get_int3()
              << std::endl;
    std::cout << "Received ints: " << recv_custom_object.get_int1() << " "
              << recv_custom_object.get_int2() << " "
              << recv_custom_object.get_int3() << std::endl;

    FAIL_IF_NOT(expected_custom.get_int1() == recv_custom_object.get_int1());
    FAIL_IF_NOT(expected_custom.get_int2() == recv_custom_object.get_int2());
    FAIL_IF_NOT(expected_custom.get_int3() == recv_custom_object.get_int3());

    std::cout << "Expected double 1: " << expected_custom.get_double1() << " ";
    std::cout << "Received double 1: " << recv_custom_object.get_double1();
    std::cout << std::endl;

    FAIL_IF_NOT(soft_equiv(expected_custom.get_double1(),
                           recv_custom_object.get_double1()));
    FAIL_IF_NOT(soft_equiv(expected_custom.get_double2(),
                           recv_custom_object.get_double2()));
    FAIL_IF_NOT(expected_custom.get_long1() == recv_custom_object.get_long1());
    FAIL_IF_NOT(expected_custom.get_long2() == recv_custom_object.get_long2());

  } catch (rtt_dsxx::assertion const & /*error*/) {
#ifdef C4_SCALAR
    PASSMSG("Successfully caught a ds++ exception while trying to use "
            "send_is_custom() in a C4_SCALAR build.");
#else
    FAILMSG("Encountered a ds++ exception while testing send_is_custom().");
#endif
  }

  // do the send receive again with a blocking version of custom sends and
  // receives

  // create some data to send/recv
  Custom my_custom_object_block(rtt_c4::node());

#ifdef C4_SCALAR
  // in saclar mode make the receive object the same as the on node object
  Custom recv_custom_object_block(rtt_c4::node());
#else
  // otherwise make it an invalid object
  Custom recv_custom_object_block(-1);
#endif

  if (rtt_c4::nodes() > 1) {
    // send data using blocking synchronous send. Custom sends check to make
    // sure that the type, T is the same size as its MPI type. Odd ranks send
    // first
    int recv_size = -1;
    if (rtt_c4::node() % 2) {
      rtt_c4::send_custom(&my_custom_object_block, 1, right, Custom::mpi_tag);
    } else {
      recv_size = rtt_c4::receive_custom(&recv_custom_object_block, 1, left,
                                         Custom::mpi_tag);
    }
    if (!(rtt_c4::node() % 2)) {
      rtt_c4::send_custom(&my_custom_object_block, 1, right, Custom::mpi_tag);

    } else {
      recv_size = rtt_c4::receive_custom(&recv_custom_object_block, 1, left,
                                         Custom::mpi_tag);
    }

    // make sure only one object was received
    if (recv_size != 1)
      ITFAILS;

    // check that the exected results match the custom type from the left rank
    Custom expected_custom(left);

    std::cout << "Expected ints: " << expected_custom.get_int1() << " "
              << expected_custom.get_int2() << " " << expected_custom.get_int3()
              << std::endl;
    std::cout << "Received ints: " << recv_custom_object_block.get_int1() << " "
              << recv_custom_object_block.get_int2() << " "
              << recv_custom_object_block.get_int3() << std::endl;

    FAIL_IF_NOT(expected_custom.get_int1() ==
                recv_custom_object_block.get_int1());
    FAIL_IF_NOT(expected_custom.get_int2() ==
                recv_custom_object_block.get_int2());
    FAIL_IF_NOT(expected_custom.get_int3() ==
                recv_custom_object_block.get_int3());

    std::cout << "Expected double 1: " << expected_custom.get_double1() << " ";
    std::cout << "Received double 1: "
              << recv_custom_object_block.get_double1();
    std::cout << std::endl;

    FAIL_IF_NOT(soft_equiv(expected_custom.get_double1(),
                           recv_custom_object_block.get_double1()));
    FAIL_IF_NOT(soft_equiv(expected_custom.get_double2(),
                           recv_custom_object_block.get_double2()));
    FAIL_IF_NOT(expected_custom.get_long1() ==
                recv_custom_object_block.get_long1());
    FAIL_IF_NOT(expected_custom.get_long2() ==
                recv_custom_object_block.get_long2());
  }
  return;
}

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_simple(ut);
    test_send_custom(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// end of tstsend_is.cc
//----------------------------------------------------------------------------//
