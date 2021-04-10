//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/test/tstquick_index.cc
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  KDE function tests
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "kde/quick_index.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <numeric>

using namespace rtt_dsxx;
using namespace rtt_c4;
using namespace rtt_kde;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
//
void test_replication(ParallelUnitTest &ut) {
  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }
    // in rep mode the max window size does nothing so set it large
    const double max_window_size = 100.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = false;
    quick_index<1> qindex = quick_index<1>(position_array, max_window_size, bins_per_dim, dd);
    // Check public data
    //  if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
    //    ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("quick_index checks pass");
  } else {
    FAILMSG("quick_index checks failed");
  }
}

void test_decomposition(ParallelUnitTest &ut) {
  if (rtt_c4::nodes() != 3)
    ITFAILS;

  int local_size = 3;
  // give the odd size to the final rank to make striding easy
  if (rtt_c4::node() == 2)
    local_size = 4;

  {
    std::vector<double> data{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::vector<std::array<double, 3>> position_array(10, std::array<double, 3>{0.0, 0.0, 0.0});
    for (int i = 0; i < 10; i++) {
      position_array[i][0] = i < 5 ? i % 5 : i % 5 + 0.5;
      position_array[i][1] = i < 5 ? 0.5 : -0.5;
    }

    // map to dd arrays with simple stride
    std::vector<double> dd_data(local_size, 0.0);
    std::vector<std::array<double, 3>> dd_position_array(local_size,
                                                         std::array<double, 3>{0.0, 0.0, 0.0});

    for (int i = 0; i < local_size; i++) {
      dd_data[i] = data[i + rtt_c4::node() * 3];
      dd_position_array[i] = position_array[i + rtt_c4::node() * 3];
    }

    // in dd mode the max window size determines the number of ghost cells
    const double max_window_size = 1.0;
    const size_t bins_per_dim = 10UL;
    const bool dd = true;
    quick_index<1> qindex = quick_index<1>(dd_position_array, max_window_size, bins_per_dim, dd);
    // Check the reuslts
    //  if (!rtt_dsxx::soft_equiv(smooth_result[i], 0.1))
    //     ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("quick_index DD checks pass");
  } else {
    FAILMSG("quick_index DD checks failed");
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ParallelUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test_replication(ut);
    if (nodes() == 3)
      test_decomposition(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstquick_index.cc
//------------------------------------------------------------------------------------------------//
