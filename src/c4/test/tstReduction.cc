//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstReduction.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 15:41:00 2002
 * \brief  C4 Reduction test.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;

using rtt_c4::C4_Req;
using rtt_c4::global_isum;
using rtt_c4::global_max;
using rtt_c4::global_min;
using rtt_c4::global_prod;
using rtt_c4::global_sum;
using rtt_c4::prefix_sum;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void elemental_reduction(rtt_dsxx::UnitTest &ut) {
  // test ints with blocking and non-blocking sums
  int xint = rtt_c4::node() + 1;
  global_sum(xint);

  int xint_send = rtt_c4::node() + 1;
  int xint_recv = 0;
  C4_Req int_request;
  global_isum(xint_send, xint_recv, int_request);
  int_request.wait();

  int int_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    int_answer += i + 1;

  if (xint != int_answer)
    ITFAILS;

  if (xint_recv != int_answer)
    ITFAILS;
  if (!rtt_c4::node())
    cout << "int: Global non-blocking sum: " << xint_recv
         << " answer: " << int_answer << endl;

  // Test with deprecated form of global_sum
  xint = rtt_c4::node() + 1;
  global_sum(xint);
  if (xint != int_answer)
    ITFAILS;

  // test longs for blocking and non-blocking sums
  long const max_long(std::numeric_limits<long>::max());
  int64_t const ten_billion(10000000000L); // 1e10 > MAX_INT
  int32_t const one_billion(1000000000L);  // 1e9 < MAX_INT

  long xlong =
      rtt_c4::node() + (max_long > ten_billion ? ten_billion : one_billion);
  global_sum(xlong);

  long xlong_send =
      rtt_c4::node() + (max_long > ten_billion ? ten_billion : one_billion);
  long xlong_recv = 0;
  C4_Req long_request;
  global_isum(xlong_send, xlong_recv, long_request);
  long_request.wait();

  long long_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    long_answer += i + (max_long > ten_billion ? ten_billion : one_billion);

  if (xlong != long_answer)
    ITFAILS;

  if (xlong_recv != long_answer)
    ITFAILS;

  // test doubles for blocking and non-blocking sums
  double xdbl = static_cast<double>(rtt_c4::node()) + 0.1;
  global_sum(xdbl);

  double xdouble_send = static_cast<double>(rtt_c4::node()) + 0.1;
  double xdouble_recv = 0;
  C4_Req double_request;
  global_isum(xdouble_send, xdouble_recv, double_request);
  double_request.wait();

  double dbl_answer = 0.0;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    dbl_answer += static_cast<double>(i) + 0.1;

  if (!soft_equiv(xdbl, dbl_answer))
    ITFAILS;

  if (!soft_equiv(xdouble_recv, dbl_answer))
    ITFAILS;

  // test product
  xlong = rtt_c4::node() + 1;
  global_prod(xlong);

  long_answer = 1;
  for (int i = 0; i < rtt_c4::nodes(); i++)
    long_answer *= (i + 1);

  if (xlong != long_answer)
    ITFAILS;

  // Test with deprecated form of global_prod
  xlong = rtt_c4::node() + 1;
  global_prod(xlong);
  if (xlong != long_answer)
    ITFAILS;

  // test min
  xdbl = 0.5 + rtt_c4::node();
  global_min(xdbl);

  if (!soft_equiv(xdbl, 0.5))
    ITFAILS;

  // Test with deprecated form of global_min
  xdbl = rtt_c4::node() + 0.5;
  global_min(xdbl);
  if (!soft_equiv(xdbl, 0.5))
    ITFAILS;

  // test max
  xdbl = 0.7 + rtt_c4::node();
  global_max(xdbl);

  if (!soft_equiv(xdbl, rtt_c4::nodes() - 0.3))
    ITFAILS;

  // Test with deprecated form of global_max
  xdbl = 0.7 + rtt_c4::node();
  global_max(xdbl);
  if (!soft_equiv(xdbl, rtt_c4::nodes() - 0.3))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Elemental reductions ok.");
  return;
}

//---------------------------------------------------------------------------//
void array_reduction(rtt_dsxx::UnitTest &ut) {
  // make a vector of doubles
  vector<double> x(100);
  vector<double> prod(100, 1.0);
  vector<double> sum(100, 0.0);
  vector<double> lmin(100, 0.0);
  vector<double> lmax(100, 0.0);

  // fill it
  for (int i = 0; i < 100; i++) {
    x[i] = rtt_c4::node() + 0.11;
    for (int j = 0; j < rtt_c4::nodes(); j++) {
      sum[i] += (j + 0.11);
      prod[i] *= (j + 0.11);
    }
    lmin[i] = 0.11;
    lmax[i] = rtt_c4::nodes() + 0.11 - 1.0;
  }

  vector<double> c;

  {
    c = x;
    global_sum(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), sum.begin(), sum.end()))
      ITFAILS;

    c = x;
    global_prod(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), prod.begin(), prod.end()))
      ITFAILS;

    c = x;
    global_min(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), lmin.begin(), lmin.end()))
      ITFAILS;

    c = x;
    global_max(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), lmax.begin(), lmax.end()))
      ITFAILS;
  }

  // Test using deprecated forms of global_sum, global_min, global_max and
  // global_prod.

  {
    c = x;
    global_sum(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), sum.begin(), sum.end()))
      ITFAILS;

    c = x;
    global_prod(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), prod.begin(), prod.end()))
      ITFAILS;

    c = x;
    global_min(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), lmin.begin(), lmin.end()))
      ITFAILS;

    c = x;
    global_max(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), lmax.begin(), lmax.end()))
      ITFAILS;
  }

  if (ut.numFails == 0)
    PASSMSG("Array reductions ok.");
  return;
}

//---------------------------------------------------------------------------//
void test_prefix_sum(rtt_dsxx::UnitTest &ut) {

  // Calculate prefix sums on rank ID with MPI call and by hand and compare the
  // output. The prefix sum on a node includes all previous node's value and the
  // value of the current node

  // test ints
  int xint = rtt_c4::node();
  int xint_prefix_sum = prefix_sum(xint);

  int int_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i < rtt_c4::node())
      int_answer += i + 1;
  }

  std::cout << "int: Prefix sum on this node: " << xint_prefix_sum;
  std::cout << " Answer: " << int_answer << std::endl;

  if (xint_prefix_sum != int_answer)
    ITFAILS;

  // test unsigned ints (start at max of signed int)
  uint32_t xuint = rtt_c4::node();
  if (rtt_c4::node() == 0)
    xuint = std::numeric_limits<int>::max();
  uint32_t xuint_prefix_sum = prefix_sum(xuint);

  uint32_t uint_answer = std::numeric_limits<int>::max();
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i < rtt_c4::node())
      uint_answer += i + 1;
  }

  std::cout << "uint32_t: Prefix sum on this node: " << xuint_prefix_sum;
  std::cout << " Answer: " << uint_answer << std::endl;

  if (xuint_prefix_sum != uint_answer)
    ITFAILS;

  // test longs
  long xlong = rtt_c4::node() + 1000;
  long xlong_prefix_sum = prefix_sum(xlong);

  long long_answer = 0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i <= rtt_c4::node() || i == 0)
      long_answer += i + 1000;
  }

  std::cout << "long: Prefix sum on this node: " << xlong_prefix_sum;
  std::cout << " Answer: " << long_answer << std::endl;

  if (xlong_prefix_sum != long_answer)
    ITFAILS;

  // test unsigned longs (start at max of unsigned int)
  uint64_t xulong = rtt_c4::node();
  if (rtt_c4::node() == 0)
    xulong = std::numeric_limits<uint32_t>::max();
  uint64_t xulong_prefix_sum = prefix_sum(xulong);

  uint64_t ulong_answer = std::numeric_limits<uint32_t>::max();
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i < rtt_c4::node())
      ulong_answer += i + 1;
  }

  std::cout << "uint64_t: Prefix sum on this node: " << xulong_prefix_sum;
  std::cout << " Answer: " << ulong_answer << std::endl;

  if (xulong_prefix_sum != ulong_answer)
    ITFAILS;

  // test floats
  float xfloat = static_cast<float>(rtt_c4::node() + 0.01);
  float xfloat_prefix_sum = prefix_sum(xfloat);

  float float_answer = 0.0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i <= rtt_c4::node() || i == 0)
      float_answer += static_cast<float>(i + 0.01);
  }

  std::cout << "float: Prefix sum on this node: " << xfloat_prefix_sum;
  std::cout << " Answer: " << float_answer << std::endl;

  if (!soft_equiv(xfloat_prefix_sum, float_answer))
    ITFAILS;

  // test doubles
  double xdbl = static_cast<double>(rtt_c4::node()) + 1.0e-9;
  double xdbl_prefix_sum = prefix_sum(xdbl);

  double dbl_answer = 0.0;
  for (int i = 0; i < rtt_c4::nodes(); i++) {
    if (i <= rtt_c4::node() || i == 0)
      dbl_answer += static_cast<double>(i) + 1.0e-9;
  }

  std::cout.precision(16);
  std::cout << "double: Prefix sum on this node: " << xdbl_prefix_sum;
  std::cout << " Answer: " << dbl_answer << std::endl;

  if (!soft_equiv(xdbl_prefix_sum, dbl_answer))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("Prefix sum ok.");
  return;
}

//---------------------------------------------------------------------------//
void test_array_prefix_sum(rtt_dsxx::UnitTest &ut) {

  // Calculate prefix sums on rank ID with MPI call and by hand and compare the
  // output. The prefix sum on a node includes all previous node's value and
  // the value of the current node

  const int array_size = 12;

  // test ints
  vector<int> xint(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xint[i] = rtt_c4::node() * 10 + i;

  prefix_sum(&xint[0], array_size);

  vector<int> int_answer(array_size, 0);
  for (int i = 0; i < array_size; ++i) {
    for (int r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        int_answer[i] += r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < xint.size(); ++i) {
    std::cout << "int: Prefix sum on this node: " << xint[i];
    std::cout << " Answer: " << int_answer[i] << std::endl;
    if (xint[i] != int_answer[i])
      ITFAILS;
  }

  // test unsigned ints (use the maximum int value to make sure all types are
  // handled correctly in the calls)
  vector<uint32_t> xuint(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xuint[i] = std::numeric_limits<int>::max() + rtt_c4::node() * 10 + i;

  prefix_sum(&xuint[0], array_size);

  vector<uint32_t> uint_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        uint_answer[i] += std::numeric_limits<int>::max() + r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < xuint.size(); ++i) {
    std::cout << "uint32_t: Prefix sum on this node: " << xuint[i];
    std::cout << " Answer: " << uint_answer[i] << std::endl;
    if (xuint[i] != uint_answer[i])
      ITFAILS;
  }

  // test long ints (use the maximum uint32_t value to make sure all types are
  // handled correctly in the calls)
  vector<int64_t> xlong(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xlong[i] = std::numeric_limits<uint32_t>::max() + rtt_c4::node() * 10 + i;

  prefix_sum(&xlong[0], array_size);

  vector<int64_t> long_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        long_answer[i] += std::numeric_limits<uint32_t>::max() + r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < xlong.size(); ++i) {
    std::cout << "int64_t: Prefix sum on this node: " << xlong[i];
    std::cout << " Answer: " << long_answer[i] << std::endl;
    if (xlong[i] != long_answer[i])
      ITFAILS;
  }

  // test unsigned long ints (use the maximum int64_t value to make sure all
  // types are handled correctly in the calls)
  vector<uint64_t> xulong(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xulong[i] = std::numeric_limits<int64_t>::max() + rtt_c4::node() * 10 + i;

  prefix_sum(&xulong[0], array_size);

  vector<uint64_t> ulong_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        ulong_answer[i] += std::numeric_limits<int64_t>::max() + r * 10 + i;
    }
  }

  for (uint32_t i = 0; i < xulong.size(); ++i) {
    std::cout << "uint64_t: Prefix sum on this node: " << xulong[i];
    std::cout << " Answer: " << ulong_answer[i] << std::endl;
    if (xulong[i] != ulong_answer[i])
      ITFAILS;
  }

  // test floats
  vector<float> xfloat(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xfloat[i] = static_cast<float>(rtt_c4::node() * 9.99 + i);

  prefix_sum(&xfloat[0], array_size);

  vector<float> float_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        float_answer[i] += static_cast<float>(r * 9.99 + i);
    }
  }

  // comparison between floats after operations needs soft_equiv with a loose
  // tolerance
  for (uint32_t i = 0; i < xfloat.size(); ++i) {
    std::cout << "float: Prefix sum on this node: " << xfloat[i];
    std::cout << " Answer: " << float_answer[i] << std::endl;
    if (!soft_equiv(xfloat[i], float_answer[i], float(1.0e-6)))
      ITFAILS;
  }

  // test doubles, try to express precision beyond float to test type handling
  // in function calls
  vector<double> xdouble(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i)
    xdouble[i] = rtt_c4::node() * 9.000000000002 + i;

  prefix_sum(&xdouble[0], array_size);

  vector<double> double_answer(array_size, 0);
  for (int32_t i = 0; i < array_size; ++i) {
    for (int32_t r = 0; r < rtt_c4::nodes(); ++r) {
      if (r <= rtt_c4::node())
        double_answer[i] += r * 9.000000000002 + i;
    }
  }

  for (uint32_t i = 0; i < xdouble.size(); ++i) {
    std::cout << "double: Prefix sum on this node: " << xdouble[i];
    std::cout << " Answer: " << double_answer[i] << std::endl;
    if (!soft_equiv(xdouble[i], double_answer[i]))
      ITFAILS;
  }

  if (ut.numFails == 0)
    PASSMSG("Array prefix sum ok.");
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    elemental_reduction(ut);
    array_reduction(ut);
    test_prefix_sum(ut);
    test_array_prefix_sum(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstReduction.cc
//---------------------------------------------------------------------------//
