//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ode/test/tstFunction_Traits.cc
 * \author Kent Budge
 * \date   Wed Aug 18 10:28:16 2004
 * \brief  
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ode/Function_Traits.hh"

#include <iostream>
#include <typeinfo>

using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

class Test_Functor {
public:
  typedef double return_type;
};

void tstFunction_Traits(UnitTest &ut) {
  if (typeid(Function_Traits<double (*)(double)>::return_type) !=
      typeid(double))
    ut.failure("return_type NOT correct");
  else
    ut.passes("return_type correct");

  if (typeid(Function_Traits<Test_Functor>::return_type) !=
      typeid(Test_Functor::return_type))
    ut.failure("return_type NOT correct");
  else
    ut.passes("return_type correct");
  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstFunction_Traits(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstFunction_Traits.cc
//---------------------------------------------------------------------------//
