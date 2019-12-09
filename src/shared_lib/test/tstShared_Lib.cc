//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstShared_Lib.cc
 * \author Rob Lowrie
 * \date   Thu Apr 15 23:03:32 2004
 * \brief  Tests Shared_Lib
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../Shared_Lib.hh"
#include "shared_lib_test.hh"
#include "ds++/Release.hh"
#include <ds++/Assert.hh>
#include <ds++/Soft_Equivalence.hh>

#include "Foo_Base.hh"

using namespace std;
using namespace rtt_shared_lib;
using rtt_shared_lib_test::Foo_Base;

using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

bool test_supported() {
  bool supported = Shared_Lib::is_supported();

  // Make sure we can't create a Shared_Lib on unsupported platforms.

  if (not supported) {
    bool caught = false;

    try {
      Shared_Lib s;
    } catch (const rtt_dsxx::assertion &ass) {
      caught = true;
      std::ostringstream m;
      m << "Excellent! Caught assertion for unsupported platforms.";
      PASSMSG(m.str());
    }

    if (not caught)
      ITFAILS;
  }

  // Done testing.

  if (rtt_shared_lib_test::passed) {
    if (supported) {
      PASSMSG("test_supported() ok, platform supported.");
    } else {
      PASSMSG("test_supported() ok, platform unsupported.");
    }
  }

  return supported;
}

void test_simple() {
  string so_name("foo.so");

  Shared_Lib so_lib(so_name);

  if (so_lib.get_file_name() != so_name)
    ITFAILS;
  if (!so_lib.is_open())
    ITFAILS;

  // Creator function pointer for Foo objects.
  typedef Foo_Base *(*Creator)(double);

  // Destroyer function pointer for Foo objects.
  typedef void (*Destroyer)(Foo_Base *);

  // Grab the creator function.
  Creator my_creator = so_lib.get_function<Creator>("my_creator");

  // Grab the destroyer function.
  Destroyer my_destroyer = so_lib.get_function<Destroyer>("my_destroyer");

  const double base = 4.1;
  const double x = 2.34342;

  // Create the Foo object and check it.
  Foo_Base *foo = my_creator(base);
  if (!foo)
    ITFAILS;
  if (!soft_equiv(foo->compute(x), x * x + base))
    ITFAILS;

  { // Test copy ctor
    Shared_Lib sc(so_lib);

    if (sc.get_file_name() != so_name)
      ITFAILS;
    Creator c = sc.get_function<Creator>("my_creator");
    Destroyer d = sc.get_function<Destroyer>("my_destroyer");
    Foo_Base *foo2 = c(base);
    if (!foo2)
      ITFAILS;
    if (!soft_equiv(foo2->compute(x), x * x + base))
      ITFAILS;
    d(foo2);
  }

  { // Test assignment
    Shared_Lib sc;

    sc = so_lib;

    if (sc.get_file_name() != so_name)
      ITFAILS;
    Creator c = sc.get_function<Creator>("my_creator");
    Destroyer d = sc.get_function<Destroyer>("my_destroyer");
    Foo_Base *foo2 = c(base);
    if (!foo2)
      ITFAILS;
    if (!soft_equiv(foo2->compute(x), x * x + base))
      ITFAILS;
    d(foo2);
  }

  { // check get_function of a function not in the library
    bool caught = false;
    try {
      /* Creator no = */ so_lib.get_function<Creator>("not_in_lib");
    } catch (const rtt_dsxx::assertion &ass) {
      caught = true;
      std::ostringstream m;
      m << "Excellent! Caught assertion for get_function().";
      PASSMSG(m.str());
    }

    if (not caught)
      ITFAILS;
  }

  // Done with so_lib, so let's destroy foo and close so_lib.
  my_destroyer(foo);
  so_lib.close();
  if (so_lib.is_open())
    ITFAILS;

#ifdef REQUIRE_ON
  { // check open() of an empty file name
    bool caught = false;
    try {
      Shared_Lib t;
      t.open("");
    } catch (const rtt_dsxx::assertion &ass) {
      caught = true;
      std::ostringstream m;
      m << "Excellent! Caught assertion for open().";
      PASSMSG(m.str());
    }

    if (not caught)
      ITFAILS;
  }
#endif

  // Done testing

  if (rtt_shared_lib_test::passed) {
    PASSMSG("test_simple() ok.");
  }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  // version tag
  for (int arg = 1; arg < argc; arg++)
    if (std::string(argv[arg]) == "--version") {
      std::cout << argv[0] << ": version " << rtt_dsxx::release() << std::endl;
      return 0;
    }

  try {
    // >>> UNIT TESTS

    if (test_supported()) {
      test_simple();
    }
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing tstShared_Lib, " << err.what()
              << std::endl;
    return 1;
  } catch (...) {
    std::cout << "ERROR: While testing tstShared_Lib, "
              << "An unknown exception was thrown." << std::endl;
    return 1;
  }

  // status of test
  std::cout << std::endl;
  std::cout << "*********************************************" << std::endl;
  if (rtt_shared_lib_test::passed) {
    std::cout << "**** tstShared_Lib Test: PASSED" << std::endl;
  }
  std::cout << "*********************************************" << std::endl;
  std::cout << std::endl;

  std::cout << "Done testing tstShared_Lib." << std::endl;
  return 0;
}

//---------------------------------------------------------------------------//
// end of tstShared_Lib.cc
//---------------------------------------------------------------------------//
