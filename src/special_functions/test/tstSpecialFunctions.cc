//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/test/test_sf.cc
 * \author Kelly Thompson
 * \date   Tue Sep 27 12:49:39 2005
 * \brief  Unit tests for kronecker_delta and factorial.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "special_functions/Factorial.hh"
#include "special_functions/KroneckerDelta.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstKdelta(rtt_dsxx::UnitTest &ut) {
  using rtt_sf::kronecker_delta;
  if (kronecker_delta(0, 0) == 1)
    ut.passes("Found kronecker_delta(0,0) == 1, kronecker_delta is working.");
  else
    ut.failure(
        "Found kronecker_delta(0,0) != 1, kronecker_delta is not working.");
  if (kronecker_delta(0, 1) == 0)
    ut.passes("Found kronecker_delta(0,1) == 0, kronecker_delta is working.");
  else
    ut.failure(
        "Found kronecker_delta(0,1) != 0, kronecker_delta is not working.");
  if (kronecker_delta(1, 1) == 1)
    ut.passes("Found kronecker_delta(1,1) == 1, kronecker_delta is working.");
  else
    ut.failure(
        "Found kronecker_delta(1,1) != 1, kronecker_delta is not working.");
  if (kronecker_delta(1, 0) == 0)
    ut.passes("Found kronecker_delta(1,0) == 0, kronecker_delta is working.");
  else
    ut.failure(
        "Found kronecker_delta(1,0) != 0, kronecker_delta is not working.");
  if (kronecker_delta(-1, 0) == 0)
    ut.passes("Found kronecker_delta(-1,0) == 0, kronecker_delta is working.");
  else
    ut.failure(
        "Found kronecker_delta(-1,0) != 0, kronecker_delta is not working.");
  if (kronecker_delta(-1, -1) == 1)
    ut.passes("Found kronecker_delta(-1,-1) == 1, kronecker_delta is working.");
  else
    ut.failure(
        "Found kronecker_delta(-1,-1) != 1, kronecker_delta is not working.");

  unsigned uZero(0);
  unsigned uOne(1);
  long lZero(0);
  long lOne(1);

  if (kronecker_delta(uOne, uZero) == uZero)
    ut.passes("Found kronecker_delta<unsigned>(uOne,uZero) == uZero, "
              "kronecker_delta is working.");
  else
    ut.failure("Found kronecker_delta<unsigned>(uOne,uZero) != uZero, "
               "kronecker_delta is not working.");

  if (kronecker_delta(lOne, lZero) == lZero)
    ut.passes("Found kronecker_delta<long>(uOne,uZero) == uZero, "
              "kronecker_delta is working.");
  else
    ut.failure("Found kronecker_delta<long>(uOne,uZero) != uZero, "
               "kronecker_delta is not working.");

  return;
}

//---------------------------------------------------------------------------//

void tstFactorial(rtt_dsxx::UnitTest &ut) {
  using rtt_sf::factorial;

  // Test factorial

  if (factorial(0) == 1)
    ut.passes("Found factorial(0) == 1, factorial is working.");
  else
    ut.failure("Found factorial(0) != 1, factorial is not working.");
  if (factorial(1) == 1)
    ut.passes("Found factorial(1) == 1, factorial is working.");
  else
    ut.failure("Found factorial(1) != 1, factorial is not working.");
  if (factorial(2) == 2)
    ut.passes("Found factorial(2) == 2, factorial is working.");
  else
    ut.failure("Found factorial(2) != 2, factorial is not working.");
  if (factorial(3) == 6)
    ut.passes("Found factorial(3) == 6, factorial is working.");
  else
    ut.failure("Found factorial(3) != 6, factorial is not working.");

  try {
    factorial(13);
    ut.failure("factorial(13) failed to trigger out of range error");
  } catch (std::range_error &) {
    ut.passes("factorial(13) correctly triggered out of range error");
  }

  if (factorial(-3) == 1)
    ut.passes("Found factorial(-3) == 1, factorial is working.");
  else
    ut.failure("Found factorial(-3) != 1, factorial is not working.");

  unsigned uOne(1);
  long lOne(1);

  if (factorial(uOne) == uOne)
    ut.passes(
        "Found factorial<unsigned>(1) == unsigned(1), factorial is working.");
  else
    ut.failure("Found factorial<unsigned>(1) != unsigned(1), factorial is not "
               "working.");
  if (factorial(lOne) == lOne)
    ut.passes("Found factorial<long>(1) == long(1), factorial is working.");
  else
    ut.failure(
        "Found factorial<long>(1) != long(1), factorial is not working.");
  return;
}
//---------------------------------------------------------------------------//

void tstFF(rtt_dsxx::UnitTest &ut) {
  using rtt_sf::factorial_fraction;

  // Test factorial_fraction

  if (factorial_fraction(1, 1) == 1)
    ut.passes("factorial_fraction(1,1) == 1");
  else
    ut.failure("factorial_fraction(1,1) != 1");

  if (factorial_fraction(6, 4) == 30)
    ut.passes("factorial_fraction(6,4) == 30");
  else
    ut.failure("factorial_fraction(6,4) != 30");

  if (rtt_dsxx::soft_equiv(factorial_fraction(1, 2), 0.5))
    ut.passes("factorial_fraction(1,2) == 0.5");
  else
    ut.failure("factorial_fraction(1,2) != 0.5");

  if (factorial_fraction(100, 99) == 100)
    ut.passes("factorial_fraction(100,99) == 100");
  else
    ut.failure("factorial_fraction(100,99) != 100");

  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  using namespace rtt_sf;
  using namespace std;
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tstFactorial(ut);
    tstKdelta(ut);
    tstFF(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of test_sf.cc
//---------------------------------------------------------------------------//
