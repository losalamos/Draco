//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FortranCheck/test/fc_derived_type.cc
 * \author Allan Wollaber
 * \date   Tue Jul 10 12:48:13 MDT 2012
 * \brief  Test Fortran main calling C with a derived type
 * \note   Copyright (c) 2012 Los Alamos National Security, LLC
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//
#include <stdint.h>
#include <cmath>
#include <iostream>

/// Define the derived type as a C structure
struct my_informative_type
{
   double  some_double;
   int     some_int;
   int64_t some_large_int;
};


// A simple function to test for valid values in a Fortran derived type
extern "C" 
void rtt_test_derived_type(const my_informative_type& mit, int& error_code)
{
    std::cout << "On C-interface, derived type has double = " <<
                  mit.some_double << ", int = " << mit.some_int
                  << ", large_int = " << mit.some_large_int << std::endl;
    std::cout << std::endl;

    error_code = 0;

    if (std::abs(mit.some_double - 3.141592654) > 1e-9)
    {
       error_code = 1;
       return;
    }
    else if (mit.some_int != 137)
    {
       error_code = 2;
       return;
    }
    else if (mit.some_large_int != ( (2LL) << 33) )
       error_code = 3;

    return;
}
