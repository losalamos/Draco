//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/tstEnsight_Stream.cc
 * \author Rob Lowrie
 * \date   Fri Nov 12 22:52:46 2004
 * \brief  Test for Ensight_Stream.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Packing_Utils.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "viz/Ensight_Stream.hh"

using namespace std;
using rtt_viz::Ensight_Stream;

//---------------------------------------------------------------------------//
// Utility functions
//---------------------------------------------------------------------------//

// Reads binary value from stream.
template <typename T> void binary_read(ifstream &stream, T &v) {
  char *vc = new char[sizeof(T)];
  stream.read(vc, sizeof(T));

  rtt_dsxx::Unpacker p;
  p.set_buffer(sizeof(T), vc);
  p.unpack(v);

  delete[] vc;
}

// Various overloaded read functions.

void readit(ifstream &stream, const bool binary, double &d) {
  if (binary) {
    float x;
    binary_read(stream, x);
    d = x;
  } else
    stream >> d;
}

void readit(ifstream &stream, const bool binary, int &d) {
  if (binary)
    binary_read(stream, d);
  else
    stream >> d;
}

void readit(ifstream &stream, const bool binary, string &s) {
  if (binary) {
    s.resize(80);
    for (int i = 0; i < 80; ++i)
      stream.read(&s[i], 1);
  } else
    stream >> s;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_simple(rtt_dsxx::UnitTest &ut, bool const binary) {
  // Dump a few values into the stream

  const int i(20323);
  const string s("dog");
  const double d(112.3);
  const string file("ensight_stream.out");

  {
    Ensight_Stream f(file, binary);

    f << i << rtt_viz::endl;
    f << d << rtt_viz::endl;
    f << s << rtt_viz::endl;
  }

  // Read the file back in and check the values.

  std::ios::openmode mode = std::ios::in;

  if (binary) {
    cout << "Testing binary mode." << endl;
    mode = mode | std::ios::binary;
  } else
    cout << "Testing ascii mode." << endl;

  ifstream in(file.c_str(), mode);

  int i_in;
  readit(in, binary, i_in);
  if (i != i_in)
    ITFAILS;

  double d_in;
  readit(in, binary, d_in);
  // floats are inaccurate
  if (!rtt_dsxx::soft_equiv(d, d_in, 0.01))
    ITFAILS;

  string s_in;
  readit(in, binary, s_in);
  for (size_t k = 0; k < s.size(); ++k)
    if (s[k] != s_in[k])
      ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("test_simple() completed successfully.");
  else
    FAILMSG("test_simple() did not complet successfully.");

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {                     // >>> UNIT TESTS
    test_simple(ut, true);  // test binary
    test_simple(ut, false); // test ascii
    test_simple(ut, true);  // test binary again
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstEnsight_Stream.cc
//---------------------------------------------------------------------------//
