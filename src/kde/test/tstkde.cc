//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/test/tstkde.cc
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  <start>
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "kde/kde.hh"

using namespace rtt_dsxx;
using namespace kde;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
//
void test_instansiation(ScalarUnitTest &ut) {
  kde::kde<kde_coordinates::CART> test_kde;
  std::vector<double> data(10, 0.0);
  test_kde.reconstruction<1>(data, data, data, data);
  test_kde.reconstruction<2>(data, data, data, data);
  test_kde.reconstruction<3>(data, data, data, data);
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test_instansiation(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstkde.cc
//------------------------------------------------------------------------------------------------//

