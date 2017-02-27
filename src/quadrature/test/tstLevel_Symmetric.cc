//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tstLevel_Symmetric.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id: template_test.cc 5830 2011-05-05 19:43:43Z kellyt $
//---------------------------------------------------------------------------//

#include "quadrature_test.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "quadrature/Level_Symmetric.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    for (unsigned N = 2; N < 18; N += 2)
    // Pathological in curvilinear geometry at N>=18
    {
      cout << "Order " << N << ':' << endl;
      Level_Symmetric quadrature(N);
      if (quadrature.sn_order() != N) {
        ut.failure("NOT correct SN order");
      }
      quadrature_test(ut, quadrature);
    }
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstLevel_Symmetric.cc
//---------------------------------------------------------------------------//
