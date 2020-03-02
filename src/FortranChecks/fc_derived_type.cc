//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   FortranChecks/fc_derived_type.cc
 * \author Allan Wollaber
 * \date   Tue Jul 10 12:48:13 MDT 2012
 * \brief  Test Fortran main calling C with a derived type
 * \note   Copyright (c) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "ds++/config.h"
#include <cmath>
#include <cstdint>
#include <iostream>

//! Define the derived type as a C structure
struct my_informative_type {
  double some_double;
  int some_int;
  int64_t some_large_int;
  int *some_pointer;
  enum MG_Select { GREY = 0, MULTIGROUP = 1, ODF = 2 };
  MG_Select some_enum;
};

//! A simple function to test for valid values in a Fortran derived type
extern "C" void rtt_test_derived_type(const my_informative_type &mit,
                                      int &error_code) {
  std::cout << "In the C-interface, derived type has double = "
            << mit.some_double << std::endl
            << "int = " << mit.some_int << std::endl
            << "large_int = " << mit.some_large_int << std::endl
            << "*some_pointer[1] = " << *(mit.some_pointer) << std::endl
            << "*some_pointer[2] = " << *(mit.some_pointer + 1) << std::endl
            << "some_enum = " << mit.some_enum << std::endl;
  std::cout << std::endl;

  error_code = 0;

  if (std::abs(mit.some_double - 3.141592654) > 1e-9) {
    error_code = 1;
    return;
  } else if (mit.some_int != 137) {
    error_code = 2;
    return;
  } else if (mit.some_large_int != ((2LL) << 33)) {
    error_code = 3;
    return;
  } else if (*(mit.some_pointer) != 2003 || *(mit.some_pointer + 1) != 2012) {
    error_code = 4;
    return;
  } else if (mit.some_enum != my_informative_type::MULTIGROUP)
    error_code = 5;

  return;
}

//----------------------------------------------------------------------------//
// end of fc_derived_type.cc
//----------------------------------------------------------------------------//
