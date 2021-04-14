//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/quick_index.hh
 * \author Mathew Cleveland
 * \brief  This class generates coarse spatial indexing to quickly access
 * near-neighbor data. This additionally provides simple interpolation schemes
 * to map data to simple structured meshes. 
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_kde_quick_index_hh
#define rtt_kde_quick_index_hh

#include "c4/global.hh"
#include <array>
#include <map>
#include <vector>

namespace rtt_kde {

//================================================================================================//
/*!
 * \class quick_index
 * \brief
 *
 * Provide a hash like index of spatial distributed data along with simple
 * mapping functions.
 * 
 */
//================================================================================================//

template <size_t dim> class quick_index {
public:
  //! Default constructors.
  quick_index(const std::vector<std::array<double, 3>> &locations, const double max_window_size,
              const size_t bins_per_dimension, const bool domain_decomposed);

  //! Collect Ghost Data
  std::vector<double> collect_ghost_data(const std::vector<double> &local_data) const;

  //! Override function of 3D array ghost data.
  std::vector<std::array<double, 3>>
  collect_ghost_data(const std::vector<std::array<double, 3>> &local_data) const;

  //! Fetch list of coarse index values bound by the window
  std::vector<size_t> window_coarse_index_list(const std::array<double, 3> &window_min,
                                               const std::array<double, 3> &window_max) const;

  // PUBLIC DATA
  // Quick index initialization data
  const bool domain_decomposed;
  const size_t coarse_bin_resolution;
  const double max_window_size;
  // keep a copy of the locations
  const std::vector<std::array<double, 3>> locations;
  const size_t n_locations;

  // Global bounds
  std::array<double, 3> bounding_box_min;
  std::array<double, 3> bounding_box_max;
  // Local Data map
  std::map<size_t, std::vector<size_t>> coarse_index_map;

  // DOMAIN DECOMPOSED DATA
  // Local bounds
  std::array<double, 3> local_bounding_box_min;
  std::array<double, 3> local_bounding_box_max;
  // Ordered list of local bins (indexes values are based on the global bin structure)
  std::vector<size_t> local_bins;
  // Size of ghost data buffer
  int local_ghost_buffer_size;
  // Map used to index into a local ghost buffer
  std::map<size_t, std::vector<size_t>> local_ghost_index_map;
  // Local ghost locations (build at construction time)
  std::vector<std::array<double, 3>> local_ghost_locations;

private:
  // PRIVATE DATA
  // Map used to write local data to other processor ghost cells
  // put_window_map[global_id] = [put_rank, ghost_proc_buffer_size, ghost_proc_put_offset]
  // array is integers to accommodate mpi data types
  std::map<size_t, std::vector<std::array<int, 3>>> put_window_map;
};

} // end namespace  rtt_kde

#include "quick_index.i.hh"

#endif // rtt_kde_quick_index_hh

//------------------------------------------------------------------------------------------------//
// end of kde/quick_index.hh
//------------------------------------------------------------------------------------------------//
