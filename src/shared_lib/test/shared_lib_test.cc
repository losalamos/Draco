//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shared_lib/test/shared_lib_test.cc
 * \author Thomas M. Evans
 * \date   Wed Apr 21 14:31:07 2004
 * \brief  shared_lib testing infrastructure.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#include "shared_lib_test.hh"
#include <iostream>

namespace rtt_shared_lib_test {

//============================================================================//
// PASS/FAILURE
//============================================================================//

bool fail(int line) {
  std::cout << "Test: failed on line " << line << std::endl;
  passed = false;
  return false;
}

//----------------------------------------------------------------------------//

bool fail(int line, char *file) {
  std::cout << "Test: failed on line " << line << " in " << file << std::endl;
  passed = false;
  return false;
}

//----------------------------------------------------------------------------//

bool pass_msg(const std::string &passmsg) {
  std::cout << "Test: passed" << std::endl;
  std::cout << "     " << passmsg << std::endl;
  return true;
}

//----------------------------------------------------------------------------//

bool fail_msg(const std::string &failmsg) {
  std::cout << "Test: failed" << std::endl;
  std::cout << "     " << failmsg << std::endl;
  passed = false;
  return false;
}

//----------------------------------------------------------------------------//

void unit_test(const bool pass, int line, char *file) {
  if (pass) {
    std::cout << "Test: passed\n";
  } else {
    fail(line, file);
  }
}

//----------------------------------------------------------------------------//
// BOOLEAN PASS FLAG
//----------------------------------------------------------------------------//

bool passed = true;

} // end namespace rtt_shared_lib_test

//----------------------------------------------------------------------------//
// end of shared_lib_test.cc
//----------------------------------------------------------------------------//
