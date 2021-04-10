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

  // temp cast corse_bin_resolution to double for interpolation
  const double crd = static_cast<double>(coarse_bin_resolution);

  // build up the local hash table of into global bins
  size_t locIndex = 0;
  for (auto locItr = locations.begin(); locItr != locations.end(); locItr++) {
    auto &loc = *locItr;
    std::array<size_t, 3> index{0UL, 0UL, 0UL};
    for (size_t d = 0; d < dim; d++) {
      index[d] = static_cast<size_t>(std::floor(crd * (loc[d] - bounding_box_min[d]) /
                                                (bounding_box_max[d] - bounding_box_min[d])));
      index[d] = std::min(index[d], coarse_bin_resolution - 1);
    }
        // build up the local index hash
        const size_t global_index = index[0] + index[1] * coarse_bin_resolution +
                                    index[2] * coarse_bin_resolution * coarse_bin_resolution;
    coarse_index_map[global_index].push_back(locIndex);
    locIndex++;
  }

  // Now we need to build up ghost location map data for domain decomposed mode
  if (domain_decomposed) {
    // temporary cast of the nodes to prevent conversion warnings
    const size_t nodes = static_cast<size_t>(rtt_c4::nodes());
    const size_t node = static_cast<size_t>(rtt_c4::node());

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

    // calculate the global index range that each processor needs to
    // accommodate the specified data window size
    std::array<size_t, 3> index_min = {0UL, 0UL, 0UL};
    std::array<size_t, 3> index_max = {0UL, 0UL, 0UL};
    for (size_t d = 0; d < dim; d++) {
      // because local bounds can extend beyond the mesh we need to force a
      // positive index if necessary
      index_min[d] = static_cast<size_t>(
          std::floor(std::max(crd * (local_bounding_box_min[d] - bounding_box_min[d]) /
                                  (bounding_box_max[d] - bounding_box_min[d]),
                              0.0)));
      index_max[d] =
          static_cast<size_t>(std::floor(crd * (local_bounding_box_max[d] - bounding_box_min[d]) /
                                         (bounding_box_max[d] - bounding_box_min[d])));
      // because local bounds can extend beyond the mesh we need to floor to
      // the max bin size
      index_max[d] = std::min(index_max[d], coarse_bin_resolution - 1);
    }

    // flatten global bins
    for (size_t i = index_min[0]; i <= index_max[0]; i++) {
      for (size_t j = index_min[1]; j <= index_max[1]; j++) {
        for (size_t k = index_min[2]; k <= index_max[2]; k++) {
          size_t bin_index =
              i + j * coarse_bin_resolution + k * coarse_bin_resolution * coarse_bin_resolution;
          local_bins.push_back(bin_index);
        }
      }
    }

    // build a global map for number of entries into the global bins on each processor
    // creates a (nbins**dim)*nranks sized array
    size_t nbins = coarse_bin_resolution;
    for (size_t d = 1; d < dim; d++)
      nbins = coarse_bin_resolution;

    std::vector<int> global_index_per_bin_per_proc(nbins * nodes, 0UL);
    for (auto mapItr = coarse_index_map.begin(); mapItr != coarse_index_map.end(); mapItr++) {
      size_t gipbpp_index = mapItr->first + nbins * node;
      // must cast to an int to accomidate mpi int types.
      global_index_per_bin_per_proc[gipbpp_index] = static_cast<int>(mapItr->second.size());
    }
    rtt_c4::global_sum(&global_index_per_bin_per_proc[0], nbins * nodes);

    std::vector<int> global_need_bins_per_proc(nbins * nodes, 0UL);
    // calculate local ghost buffer size
    local_ghost_buffer_size = 0;
    for (auto binItr = local_bins.begin(); binItr != local_bins.end(); binItr++) {
      global_need_bins_per_proc[*binItr + nbins * node] += 1;
      for (size_t proc = 0; proc < nodes; proc++) {
        if (node != proc) {
          size_t gipbpp_index = *binItr + nbins * proc;
          // build up the local ghost index map
          for (int i = 0; i < global_index_per_bin_per_proc[gipbpp_index]; i++)
            local_ghost_index_map[*binItr].push_back(local_ghost_buffer_size + i);
          // accumulate the total ghost points
          local_ghost_buffer_size += global_index_per_bin_per_proc[gipbpp_index];
        }
      }
    }
    rtt_c4::global_sum(&global_need_bins_per_proc[0], nbins * nodes);

    // Build a global list of buffer sizes
    std::vector<int> proc_ghost_buffer_size(nodes, 0);
    proc_ghost_buffer_size[node] = local_ghost_buffer_size;
    rtt_c4::global_sum(&proc_ghost_buffer_size[0], nodes);

    // calculate the put map so each node knows which processor to send data
    // and where to index that data
    for (int rec_proc = 0; rec_proc < rtt_c4::nodes(); rec_proc++) {
      // calculating the offset SUCKS!!! If anyone can find a better way please help.
      int offset = 0;
      for (int send_proc = 0; send_proc < rtt_c4::node(); send_proc++) {
        if (rec_proc == send_proc)
          continue;
        for (size_t bin = 0; bin < nbins; bin++) {
          if (global_need_bins_per_proc[bin + nbins * rec_proc] > 0) {
            offset += global_index_per_bin_per_proc[bin + nbins * send_proc];
          }
        }
      }
      for (auto mapItr = coarse_index_map.begin(); mapItr != coarse_index_map.end(); mapItr++) {
        if (rtt_c4::node() != rec_proc) {
          size_t gipbpp_index = mapItr->first + nbins * rec_proc;
          if (global_need_bins_per_proc[gipbpp_index] > 0) {

            // build up map data
            put_window_map[mapItr->first].push_back(
                std::array<int, 3>{rec_proc, proc_ghost_buffer_size[rec_proc], offset});
            offset += mapItr->second.size();
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
    MPI_Win_create(&position_ghost_vector[0], local_ghost_buffer_size * sizeof(double),
                   sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    // working from my local data put the ghost data on the other ranks
    for (size_t d = 0; d < dim; d++) {
      MPI_Win_fence(0, win);
      for (auto putItr = put_window_map.begin(); putItr != put_window_map.end(); putItr++) {
        const auto &index_vector = coarse_index_map[putItr->first];
        // fill up the current ghost cell data for this dimension
        std::vector<double> put_buffer(index_vector.size());
        int putIndex = 0;
        for (auto indexItr = index_vector.begin(); indexItr != index_vector.end(); indexItr++) {
          put_buffer[putIndex] = locations[*indexItr][d];
          putIndex++;
        }
        // loop over all ranks we need to send this buffer too.
        for (auto rankItr = putItr->second.begin(); rankItr != putItr->second.end(); rankItr++) {
          std::array<int, 3> &putv = *rankItr;
          int put_rank = putv[0];
          int put_rank_buffer_size = putv[1];
          int put_offset = putv[2];
          MPI_Put(&put_buffer[0], static_cast<int>(put_buffer.size()), MPI_DOUBLE, put_rank,
                  put_offset, put_rank_buffer_size, MPI_DOUBLE, win);
        }
      }
      MPI_Win_fence(0, win);
      // alright move the position buffer to the weighting
      int posIndex = 0;
      for (auto posItr = position_ghost_vector.begin(); posItr != position_ghost_vector.end();
           posItr++) {
        local_ghost_locations[posIndex][d] = *posItr;
        posIndex++;
      }
    }
    MPI_Win_free(&win);
  } // End domain decomposed data construction
}

} // end namespace  rtt_kde

#endif // rtt_kde_quick_index_i_hh

//------------------------------------------------------------------------------------------------//
// end of kde/quick_index.i.hh
//------------------------------------------------------------------------------------------------//

