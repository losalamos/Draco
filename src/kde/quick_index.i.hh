//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/quick_index.i.hh
 * \author <user>
 * \date   <date>
 * \brief  Member definitions of class quick_index
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_kde_quick_index_i_hh
#define rtt_kde_quick_index_i_hh

#include <math.h>

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief
 * 
 * quick_index constructor. This function builds up a global indexing table to
 * quickly access data that is spatial located near each other. It breaks up
 * the data into equally spaced bins in each dimension. For domain decomposed
 * data it builds a one sided communication map to place local data that is
 * need on other processors for ghost cells. The ghost cell extents is
 * determined by the max_data_window spatial size such that any cell on the
 * local domain will have access to all points that should fall into the
 * spatial window centered on any given local point.
 *
 * \tparam dim integer specifying the data dimensionality 
 * \param[in] data locations
 * \param[in] domain_decomposed
 */
template <size_t dim>
quick_index<dim>::quick_index(const std::vector<std::array<double, 3>> &locations_,
                              const double max_window_size_, const size_t bins_per_dimension_,
                              const bool domain_decomposed_)
    : domain_decomposed(domain_decomposed_), coarse_bin_resolution(bins_per_dimension_),
      max_window_size(max_window_size_), locations(locations_), n_locations(locations_.size()) {

  // Build local bounding box
  bounding_box_min = {1e20, 1e20, 1e20};
  bounding_box_max = {-1e20, -1e20, -1e20};
  for (auto locItr = locations.begin(); locItr != locations.end(); locItr++) {
    auto &loc = *locItr;
    for (size_t d = 0; d < dim; d++) {
      if (loc[d] < bounding_box_min[d])
        bounding_box_min[d] = loc[d];
      if (loc[d] > bounding_box_max[d])
        bounding_box_max[d] = loc[d];
    }
  }

  if (domain_decomposed) {
    // Store the local bounding box and extend to maximum non-local data size
    local_bounding_box_min = bounding_box_min;
    local_bounding_box_max = bounding_box_max;
    for (size_t d = 0; d < dim; d++) {
      local_bounding_box_min[d] -= max_window_size * 0.5;
      local_bounding_box_max[d] += max_window_size * 0.5;
    }
    // Global reduce to get the global min and max
    rtt_c4::global_min(&bounding_box_min[0], 3);
    rtt_c4::global_max(&bounding_box_max[0], 3);
  }

  // build up the local hash table of into global bins
  size_t locIndex = 0;
  for (auto locItr = locations.begin(); locItr != locations.end(); locItr++) {
    auto &loc = *locItr;
    std::array<size_t, 3> index{0UL, 0UL, 0UL};
    for (size_t d = 0; d < dim; d++)
      index[d] = floor(coarse_bin_resolution * (loc[d] - bounding_box_min[d]) /
                       (bounding_box_max[d] - bounding_box_min[d]));
    // build up the local index hash
    const size_t global_index = index[0] + index[1] * coarse_bin_resolution +
                                index[2] * coarse_bin_resolution * coarse_bin_resolution;
    coarse_index_map[global_index].push_back(locIndex);
    locIndex++;
  }

  // Now we need to build up ghost location map data for domain decomposed mode
  if (domain_decomposed) {

    // calculate the global index range that each processor needs to
    // accommodate the specified data window size
    std::array<size_t, 3> index_min = {0UL, 0UL, 0UL};
    std::array<size_t, 3> index_max = {0UL, 0UL, 0UL};
    for (size_t d = 0; d < dim; d++) {
      index_min[d] =
          std::floor(coarse_bin_resolution * (local_bounding_box_min[d] - bounding_box_min[d]) /
                     (bounding_box_max[d] - bounding_box_min[d]));
      index_min[d] = std::max(index_min[d], 0UL);
      index_max[d] =
          std::floor(coarse_bin_resolution * (local_bounding_box_max[d] - bounding_box_min[d]) /
                     (bounding_box_max[d] - bounding_box_min[d]));
      index_max[d] = std::min(index_max[d], coarse_bin_resolution - 1);
    }

    // flatten global bins
    for (size_t i = index_min[0]; i <= index_max[0]; i++) {
      for (size_t j = index_min[1]; i <= index_max[1]; i++) {
        for (size_t k = index_min[2]; i <= index_max[2]; i++) {
          size_t bin_index =
              i + j * coarse_bin_resolution + k * coarse_bin_resolution * coarse_bin_resolution;
          local_bins.push_back(bin_index);
        }
      }
    }

    // build a global map for number of entries into the global bins on each processor
    // creates a (nbins**dim)*nranks sized array
    size_t nbins = pow(coarse_bin_resolution, dim);
    std::vector<size_t> global_index_per_bin_per_proc(nbins * rtt_c4::nodes(), 0UL);
    for (auto mapItr = coarse_index_map.begin(); mapItr != coarse_index_map.end(); mapItr++) {
      size_t gipbpp_index = mapItr->first + nbins * rtt_c4::node();
      global_index_per_bin_per_proc[gipbpp_index] = mapItr->second.size();
    }
    rtt_c4::global_sum(&global_index_per_bin_per_proc[0], nbins * rtt_c4::nodes());

    // calculate local ghost buffer size
    local_ghost_buffer_size = 0;
    for (auto binItr = local_bins.begin(); binItr != local_bins.end(); binItr++) {
      for (size_t proc = 0; proc < rtt_c4::nodes(); proc++) {
        if (rtt_c4::node() != proc) {
          size_t gipbpp_index = *binItr + nbins * proc;
          // build up the local ghost index map
          for (int i = 0; i < global_index_per_bin_per_proc[gipbpp_index]; i++)
            local_ghost_index_map[*binItr].push_back(local_ghost_buffer_size + i);
          // accumulate the total ghost points
          local_ghost_buffer_size += global_index_per_bin_per_proc[gipbpp_index];
        }
      }
    }

    // Build a global list of buffer sizes
    std::vector<size_t> proc_ghost_buffer_size = std::vector<size_t>(rtt_c4::nodes(), 0);
    proc_ghost_buffer_size[rtt_c4::node()] = local_ghost_buffer_size;
    rtt_c4::global_sum(&proc_ghost_buffer_size[0], rtt_c4::nodes());

    // calculate the put map so each node knows which processor to send data
    // and where to index that data
    for (auto mapItr = coarse_index_map.begin(); mapItr != coarse_index_map.end(); mapItr++) {
      for (size_t proc = 0; proc < rtt_c4::nodes(); proc++) {
        if (rtt_c4::node() != proc) {
          size_t gipbpp_index = mapItr->first + nbins * proc;
          size_t rank_count = global_index_per_bin_per_proc[gipbpp_index];
          if (rank_count > 0) {
            // calculate offset for flat_index_id vector
            size_t offset = 0;
            for (size_t sum_i = nbins * proc; sum_i < gipbpp_index; sum_i++)
              offset += global_index_per_bin_per_proc[sum_i];
            // build up map data
            put_window_map[mapItr->first] =
                std::array<size_t, 3>{proc, proc_ghost_buffer_size[proc], offset};
          }
        }
      }
    }

    // allocate the ghost data location array
    local_ghost_locations =
        std::vector<std::array<double, 3>>(local_ghost_buffer_size, {0.0, 0.0, 0.0});
    // Use one sided MPI Put commands to fill up the ghost cell location data
    std::vector<double> position_ghost_vector(local_ghost_buffer_size);
    MPI_Win win;
    MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_attach(win, &position_ghost_vector[0], local_ghost_buffer_size * sizeof(double));

    // working from my local data put the ghost data on the other ranks
    for (size_t d = 0; d < dim; d++) {
      for (auto putItr = local_ghost_index_map.begin(); putItr != local_ghost_index_map.end();
           putItr++) {
        const auto &index_vector = coarse_index_map[putItr->first];
        // fill up the current ghost cell data for this dimension
        std::vector<double> put_buffer(index_vector.size());
        int putIndex = 0;
        for (auto indexItr = index_vector.begin(); indexItr != index_vector.end(); indexItr++) {
          put_buffer[putIndex] = locations[*indexItr][d];
          putIndex++;
        }
        int put_rank = putItr->second[0];
        int put_rank_buffer_size = putItr->second[1];
        int put_offset = putItr->second[2];
        MPI_Put(&put_buffer[0], put_buffer.size(), rtt_c4::C4_Traits<double>::tag, put_rank,
                put_offset, put_rank_buffer_size, rtt_c4::C4_Traits<double>::tag, win);
      }
      // alright move the position buffer to the weighting
      rtt_c4::global_barrier(); // be sure everyone is done placing their data
      int posIndex = 0;
      for (auto posItr = position_ghost_vector.begin(); posItr != position_ghost_vector.end();
           posItr++)
        local_ghost_locations[posIndex][d] = *posItr;
    }
  } // End domain decomposed data construction
}

} // end namespace  rtt_kde

#endif // rtt_kde_quick_index_i_hh

//------------------------------------------------------------------------------------------------//
// end of kde/quick_index.i.hh
//------------------------------------------------------------------------------------------------//

