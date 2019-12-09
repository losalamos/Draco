//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstofpstream.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::determinate_swap and c4::indeterminate_swap functions
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/ofpstream.hh"
#include "ds++/Release.hh"
#include <sstream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstofpstream(UnitTest &ut) {

  unsigned const pid = rtt_c4::node();
  ofpstream out("tstofpstream.txt");

  out << "MPI rank " << pid << " reporting ..." << endl;
  out.send();
  out.shrink_to_fit();

  out << "MPI rank " << pid << " reporting a second time ..." << endl;
  out.shrink_to_fit();
  out.send();

  PASSMSG("completed serialized write without hanging or segfaulting");
}

//---------------------------------------------------------------------------//
void tstofpstream_bin(UnitTest &ut) {

  int pid = rtt_c4::node();

  // Binary write rank ids to file using ofpstream:
  {
    ofpstream out("tstofpstream.bin", std::ofstream::binary);
    out.write(reinterpret_cast<const char *>(&pid), sizeof(int));
    out.send();
    out.shrink_to_fit();
  }

  // Read file on head rank, check for correct conversion and ordering
  if (pid == 0) {
    ifstream in("tstofpstream.bin", std::ifstream::binary);
    int this_pid(-42);
    for (int a = 0; a < rtt_c4::nodes(); a++) {
      in.read(reinterpret_cast<char *>(&this_pid), sizeof(int));
      if (this_pid != a) {
        std::ostringstream msg;
        msg << "Unexpected value for this_pid = " << this_pid
            << ". Expected value a = " << a;
        FAILMSG(msg.str());
      }
    }
  }

  PASSMSG("completed serialized binary write without hanging or segfaulting");
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstofpstream(ut);
    tstofpstream_bin(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstofpstream.cc
//---------------------------------------------------------------------------//
