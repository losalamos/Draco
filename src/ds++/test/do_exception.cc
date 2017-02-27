//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/test/do_exception.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 14:33:59 2005
 * \brief  Does a floating-point exception.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Assert.hh"
#include "ds++/StackTrace.hh"
#include "ds++/fpe_trap.hh"
#include <cfenv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

//---------------------------------------------------------------------------//
/*
  Usage: do_exception test

  If test is 0, then simple floating point operations are done which should not
  cause an error.

  Otherwise, other test values should cause an exception.  Specifically, valid
  test values are

  Otherwise, other test values should cause an exception.
  Specifically, valid test values are
     1: test double division by zero
     2: test sqrt(-1)
     3: test overflow

  The file output.dat documents what happened during all tests.
*/
void run_test(int /*argc*/, char *argv[]) {

  bool const abortWithInsist(true);
  rtt_dsxx::fpe_trap fpet(abortWithInsist);
  if (fpet.enable()) {
    // Platform supported.
    cout << "\n- fpe_trap: This platform is supported.\n";
    if (!fpet.active())
      cout << "- fpe_trap: active flag set to false was not expected.\n";
  } else {
    // Platform not supported.
    cout << "\n- fpe_trap: This platform is not supported\n";
    if (fpet.active())
      cout << "- fpe_trap: active flag set to true was not expected.\n";
    return;
  }

  // Accept a command line argument with value 0, 1, 2 or 3.
  int test(-101);
  sscanf(argv[1], "%d", &test);
  Insist(test >= 0 && test <= 3, "Bad test value.");

  double zero(0.0); // for double division by zero
  double neg(-1.0); // for sqrt(-1.0)
  double result(-1.0);

  // Certain tests may be optimized away by the compiler, by recogonizing the
  // constants set above and precomputing the results below.  So do something
  // here to hopefully avoid this.  This tricks the optimizer, at least for gnu
  // and KCC.

  if (test < -100) { // this should never happen
    Insist(0, "Something is very wrong.");
    zero = neg = 1.0; // trick the optimizer?
  }

  std::ostringstream msg;
  switch (test) {
  case 0:
    // Should not throw a IEEE floating point exception.
    cout << "- Case zero: this operation should not throw a SIGFPE."
         << " The result should be 2...\n";
    result = 1.0 + zero + sqrt(-neg);
    cout << "  result = " << result << endl;
    break;
  case 1:
    cout << "- Case one: trying a div_by_zero operation...\n";
    result = 1.0 / zero; // should fail here
    cout << "  result = " << 1.0 * result;
    break;
  case 2:
    // http://en.cppreference.com/w/cpp/numeric/math/sqrt
    cout << "- Case two: trying to evaluate sqrt(-1.0)...\n";
    result = std::sqrt(neg); // should fail here
    cout << "  result = " << result;
    break;
  case 3:
    cout << "- Case three: trying to cause an overflow condition...\n";
    result = 2.0;
    std::vector<double> data;
    for (size_t i = 0; i < 100; i++) {
      // should fail at some i
      result = result * result * result * result * result;
      data.push_back(result); // force optimizer to evaluate the above line.
    }
    cout << "  result = " << result << endl;
    break;
  }
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  Insist(argc == 2, "Wrong number of args.");

  try {
    run_test(argc, argv);
  } catch (exception &err) {
    if (rtt_dsxx::fpe_trap().enable()) {
      cout << "While running " << argv[0] << ", "
           << "a SIGFPE was successfully caught.\n\t"
           << "what = " << err.what() << endl;
      return 0;
    } else {
      cout << "ERROR: While running " << argv[0] << ", "
           << "An exception was caught when it was not expected.\n\t"
           << "what = " << err.what() << endl;
    }
  } catch (...) {
    cout << "ERROR: While testing " << argv[0] << ", "
         << "An unknown exception was thrown." << endl;
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------------------//
// end of do_exception.cc
//---------------------------------------------------------------------------//
