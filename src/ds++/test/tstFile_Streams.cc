//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstFile_Streams.cc
 * \author Rob Lowrie
 * \date   Sun Nov 21 19:36:12 2004
 * \brief  Tests File_Input and File_Output.
 * \note   Copyright 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/File_Streams.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using rtt_dsxx::File_Input;
using rtt_dsxx::File_Output;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_fileio(rtt_dsxx::UnitTest &ut, const bool binary) {
  string filename("file_streams.");

  if (binary)
    filename += "binary";
  else
    filename += "ascii";

  int i = 5;
  string s = "  a string with spaces  ";
  double x = 5.6;
  bool bf = false;
  bool bt = true;

  // write the data

  {
    File_Output f(filename, binary);
    f << i;

    // here's how you write strings:
    int ssize = s.size();
    f << ssize;
    for (int k = 0; k < ssize; k++)
      f << s[k];

    f << x << bf << bt;
  }

  // read the data and make sure it's the same

  {
    int i_in;
    double x_in;
    string s_in;
    bool bf_in;
    bool bt_in;

    File_Input f(filename);
    f >> i_in;

    if (i != i_in)
      ITFAILS;

    // here's how you read strings:
    int ssize;
    f >> ssize;
    if (ssize != int(s.size()))
      ITFAILS;
    s_in.resize(ssize);
    for (int k = 0; k < ssize; k++)
      f >> s_in[k];

    if (s != s_in)
      ITFAILS;

    f >> x_in >> bf_in >> bt_in;

    if (!soft_equiv(x, x_in))
      ITFAILS;
    if (bf != bf_in)
      ITFAILS;
    if (bt != bt_in)
      ITFAILS;

    File_Input fnull("");
  }

  // test some corner cases

  {
    File_Output f;
    f.close();

    f.open("File_Stream_last_was_char.txt");
    f << 'c';
    f.close();

    File_Input fr("File_Stream_last_was_char.txt");
    char c;
    fr >> c;
    if (c != 'c')
      ITFAILS;

    fr.open("File_Stream_last_was_char.txt");
    fr >> c;
    if (c != 'c')
      ITFAILS;

    f.open("File_Stream_last_was_char.txt", false);
    f.open("File_Stream_last_was_char.txt", false);
    f << 'c';
    f.close();
    fr.open("File_Stream_last_was_char.txt");
    fr >> c;
    if (c != 'c')
      ITFAILS;
  }

  if (ut.numFails == 0) {
    ostringstream m;
    m << "test_fileio(";
    if (binary)
      m << "binary";
    else
      m << "ascii";
    m << ") ok.";
    PASSMSG(m.str());
  }
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_fileio(ut, false); // ascii
    test_fileio(ut, true);  // binary
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstFile_Streams.cc
//---------------------------------------------------------------------------//
