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
#include "ds++/dbc.hh"
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
    //------------------------
    if (qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;
    // Check global bounding box
    if (!soft_equiv(qindex.bounding_box_min[0], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], 1e20))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 1e20))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 4.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], -1e20))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], -1e20))
      ITFAILS;
    // Check local coarse_index map
    // build up a global gold to check the map
    std::map<size_t, std::vector<size_t>> gold_map;
    gold_map[0] = {0};
    gold_map[1] = {5};
    gold_map[2] = {1};
    gold_map[3] = {6};
    gold_map[4] = {2};
    gold_map[5] = {7};
    gold_map[6] = {3};
    gold_map[7] = {8};
    gold_map[8] = {4};
    gold_map[9] = {9};
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;
    for (auto mapItr = qindex.coarse_index_map.begin(); mapItr != qindex.coarse_index_map.end();
         mapItr++)
      for (size_t i = 0; i < mapItr->second.size(); i++)
        if (gold_map[mapItr->first][i] != mapItr->second[i])
          ITFAILS;
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
    // This cell spatial ordering is difficult for this setup in that every
    // rank requires a sub set of information from every other rank
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
    // Check the local state data
    if (!qindex.domain_decomposed)
      ITFAILS;
    if (qindex.coarse_bin_resolution != bins_per_dim)
      ITFAILS;
    if (!soft_equiv(qindex.max_window_size, max_window_size))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[0], 0.0))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[1], 1e20))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_min[2], 1e20))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[0], 4.5))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[1], -1e20))
      ITFAILS;
    if (!soft_equiv(qindex.bounding_box_max[2], -1e20))
      ITFAILS;
    // Check local coarse_index map
    // local indexing will not match the domain replicated case (different
    // number of points per rank so different local indexing)
    std::map<size_t, std::vector<size_t>> gold_map;
    if (rtt_c4::node() == 0) {
      gold_map[0] = {0}; // 0.0
      gold_map[2] = {1}; // 1.0
      gold_map[4] = {2}; // 2.0
    } else if (rtt_c4::node() == 1) {
      gold_map[6] = {0}; // 3.0
      gold_map[8] = {1}; // 4.0
      gold_map[1] = {2}; // 0.5
    } else {
      gold_map[3] = {0}; // 1.5
      gold_map[5] = {1}; // 2.5
      gold_map[7] = {2}; // 3.5
      gold_map[9] = {3}; // 4.5
    }
    if (gold_map.size() != qindex.coarse_index_map.size())
      ITFAILS;
    for (auto mapItr = qindex.coarse_index_map.begin(); mapItr != qindex.coarse_index_map.end();
         mapItr++)
      for (size_t i = 0; i < mapItr->second.size(); i++)
        if (gold_map[mapItr->first][i] != mapItr->second[i])
          ITFAILS;

    // Check Domain Decomposed Data
    // local bounding box extends beyond local data based on the window size
    if (rtt_c4::node() == 0) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 2.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], -1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], -1e20))
        ITFAILS;
    } else if (rtt_c4::node() == 1) {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 0.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 4.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], -1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], -1e20))
        ITFAILS;
    } else {
      if (!soft_equiv(qindex.local_bounding_box_min[0], 1.0))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[1], 1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_min[2], 1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[0], 4.5))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[1], -1e20))
        ITFAILS;
      if (!soft_equiv(qindex.local_bounding_box_max[2], -1e20))
        ITFAILS;
    }
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
