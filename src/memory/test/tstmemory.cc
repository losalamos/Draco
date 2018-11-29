//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   memory/test/tstmemory.cc
 * \author Kent G. Budge, Kelly G. Thompson
 * \brief  memory test.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "memory/memory.hh"
#include <limits>
#include <sstream>

using namespace std;
using namespace rtt_memory;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tst_memory(rtt_dsxx::UnitTest &ut) {
  if (total_allocation() == 0) {
    PASSMSG("correct initial total allocation");
  } else {
    FAILMSG("NOT correct initial total allocation");
  }
  if (peak_allocation() == 0) {
    PASSMSG("correct initial peak allocation");
  } else {
    FAILMSG("NOT correct initial peak allocation");
  }
  if (largest_allocation() == 0) {
    PASSMSG("correct initial largest allocation");
  } else {
    FAILMSG("NOT correct initial largest allocation");
  }

  set_memory_checking(true);

  double *array = new double[20];

  double *array2 = new double[30];

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == 50 * sizeof(double)) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }

  report_leaks(cerr);
#else
  PASSMSG("memory diagnostics not checked for this build");
#endif

  delete[] array;
  delete[] array2;

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == 0) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }
#endif

  // Just to try to exercise the sized delete version.
  int *scalar = new int;

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == sizeof(int)) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double) + sizeof(int)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }
#endif

  delete scalar;

#if DRACO_DIAGNOSTICS & 2
  if (total_allocation() == 0) {
    PASSMSG("correct total allocation");
  } else {
    FAILMSG("NOT correct total allocation");
  }
  if (peak_allocation() >= 50 * sizeof(double) + sizeof(int)) {
    PASSMSG("correct peak allocation");
  } else {
    FAILMSG("NOT correct peak allocation");
  }
  if (largest_allocation() >= 30 * sizeof(double)) {
    PASSMSG("correct largest allocation");
  } else {
    FAILMSG("NOT correct largest allocation");
  }
  report_leaks(cerr);
#endif
}

//---------------------------------------------------------------------------//
void tst_bad_alloc(rtt_dsxx::UnitTest &ut) {
  size_t const num_fails_start(ut.numFails);

#if DRACO_DIAGNOSTICS & 2

  std::cout << "\nTesting special handler (stack trace) for "
            << "std::bad_alloc...\n"
            << std::endl;

  // Set a specialized memory handler.
  // std::set_new_handler(rtt_memory::out_of_memory_handler);

  try {
    // trigger a std::bad_alloc exception
    std::cout << "Attempt to allocate some memory." << std::endl;
    uint64_t const new_size(static_cast<uint64_t>(100 * 1024) *
                            static_cast<uint64_t>(1024 * 1024));

    char *pBigArray = new char[new_size];

    // should never get here.
    {
      pBigArray[0] = 'a';
      std::cout << "pBigArray[0] = " << pBigArray[0] << std::endl;
      std::cout << "total_allocation = " << total_allocation() << std::endl;
    }

    delete[] pBigArray;
  } catch (std::bad_alloc &ba) {
    std::cout << "Successfully caught an expected std::bad_alloc exception."
              << std::endl;
  } catch (...) {
    FAILMSG("Failed to catch a bad_alloc.");
  }

#else // DRACO_DIAGNOSTICS & 2

#endif // DRACO_DIAGNOSTICS & 2

  if (ut.numFails > num_fails_start)
    FAILMSG("Test failures in tst_bad_alloc.");
  else
    PASSMSG("tst_bad_alloc completed sucessfully.");

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_memory(ut);
    tst_bad_alloc(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstProcmon.cc
//---------------------------------------------------------------------------//
