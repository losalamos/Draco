//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FortranCheck/test/tstDerivedTypes.cc
 * \author Kelly Thompson
 * \date   Tuesday, Jun 12, 2012, 16:03 pm
 * \brief  Test C++ main linking a Fortran library.
 * \note   Copyright (c) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// forward declaration of f90 functions
extern "C" void test_derived_types(void);

//---------------------------------------------------------------------------//
int main(int /*argc*/, char * /*argv*/ []) {
  test_derived_types();
  return 0;
}

//---------------------------------------------------------------------------//
// end of tstDerivedTypes.cc
//---------------------------------------------------------------------------//
