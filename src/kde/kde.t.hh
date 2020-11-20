//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.t.hh
 * \author Mathew Cleveland
 * \date   November 10th 2020
 * \brief  Templated KDE functions for various dimensions and coordinate
 * systems
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

// clang-format off

#ifndef kde_kde_t_hh
#define kde_kde_t_hh

#include "kde.hh"
#include <iostream>
#include <numeric>
#include <math.h>

namespace kde {

/*!
 * reconstruction 
 * \brief
 *
 * Cartesian geometry reconstruction of a 1D distribution
 *
 * \param[in] distribution 
 * \param[in] position
 * \param[in] band_width
 * \param[inout] final_distribution returned final function distribution
 *
 * Test of kde.
 */ 
template<>
template<>
std::vector<double> kde<kde_coordinates::CART>::reconstruction<1>(const
        std::vector<double> &distribution, const
        std::vector<std::array<double,3>> &position, const
        std::vector<std::array<double,3>> &band_width, const bool
        domain_decomposed) const{
    const int64_t local_size = distribution.size();
    Check(static_cast<int64_t>(position.size())==local_size);
    Check(static_cast<int64_t>(band_width.size())==local_size);
    
    int64_t size = local_size;
    int64_t global_lower_bound = 0;
    int64_t global_upper_bound = local_size;
    // minimize global values and only allocate them in DD problems
    std::vector<double> global_contribution;
    std::vector<double> global_x_position;
    double global_max_x=-1e60;
    double global_min_x=1e60;
    for (int i = 0; i<local_size; i++){
        global_max_x = std::max(global_max_x,position[i][0]);
        global_min_x = std::min(global_min_x,position[i][0]);
        Insist(!(band_width[i][0]<0.0),"bandwidth must be positive");
    }
    if(domain_decomposed){
      // calculate global off sets and min/max values
      int n_ranks = rtt_c4::nodes();
      std::vector<int64_t> rank_size(rtt_c4::nodes(),0);
      rank_size[rtt_c4::node()] = local_size;
      rtt_c4::global_sum(rank_size.data(), n_ranks);
      size = std::accumulate(rank_size.begin(),rank_size.end(),0.0);
      std::vector<int64_t> accum_rank_size(rank_size);
      std::partial_sum(rank_size.begin(), rank_size.end(),
                   accum_rank_size.begin());
      rtt_c4::global_max(global_max_x);
      rtt_c4::global_min(global_min_x);

      if(rtt_c4::node()>0){
          global_lower_bound = accum_rank_size[rtt_c4::node()-1];
          global_upper_bound = accum_rank_size[rtt_c4::node()];
      }


      // set up global arrays
      global_contribution.resize(size,0.0);
      global_x_position.resize(size,0.0);

      // build up global positions
      for (int i = 0; i<local_size; i++)
          global_x_position[i+global_lower_bound] = position[i][0];

      rtt_c4::global_sum(global_x_position.data(),size);

    }

    const double scale = std::accumulate(distribution.begin(), distribution.end(), 0.0);
    std::vector<double> result(local_size,0.0);
    std::vector<double> normal(local_size,0.0);
    std::vector<double> residual(local_size,0.0);
    std::vector<double> smooth_residual(local_size,0.0);


    // accumulate weighted contribution from all kernels for conservation correction
    for (int i = 0; i<local_size; i++){
        const double x0 = position[i][0];
        // fetch local contribution
        for (int j = 0; j<local_size; j++){
            const double x = position[j][0];
            const double h = band_width[j][0];
            const double u = (x0-x)/h;
            const double uxl = ((x-global_min_x)+(x0-global_min_x))/h;
            const double uxu = ((global_max_x-x)+(global_max_x-x0))/h;
            normal[j] += (epan_kernel(u)+epan_kernel(uxl)+epan_kernel(uxu))/h;
        }
    }
    // fetch lower ranks nonlocal contribution
    for (int i = 0; i<global_lower_bound; i++){
        const double x0 = global_x_position[i];
        // fetch local contribution
        for (int j = 0; j<local_size; j++){
            const double h = band_width[j][0];
            const double x = position[j][0];
            const double u = (x0-x)/h;
            const double uxl = ((x-global_min_x)+(x0-global_min_x))/h;
            const double uxu = ((global_max_x-x)+(global_max_x-x0))/h;
            normal[j] += (epan_kernel(u)+epan_kernel(uxl)+epan_kernel(uxu))/h;
        }
    }
    // fetch upper ranks nonlocal contribution
    for (int i = global_upper_bound; i<size; i++){
        const double x0 = global_x_position[i];
        // fetch local contribution
        for (int j = 0; j<local_size; j++){
            const double h = band_width[j][0];
            const double x = position[j][0];
            const double u = (x0-x)/h;
            const double uxl = ((x-global_min_x)+(x0-global_min_x))/h;
            const double uxu = ((global_max_x-x)+(global_max_x-x0))/h;
            normal[j] += (epan_kernel(u)+epan_kernel(uxl)+epan_kernel(uxu))/h;
        }
    }

    // now apply the kernel to the local ranks
    for (int i = 0; i<local_size; i++){
        const double x0 = position[i][0];
        // fetch local contribution
        for (int j = 0; j<local_size; j++){
            const double x = position[j][0];
            const double h = band_width[j][0];
            const double u = (x0-x)/h;
            const double uxl = ((x-global_min_x)+(x0-global_min_x))/h;
            const double uxu = ((global_max_x-x)+(global_max_x-x0))/h;
            result[i] += distribution[j]*(epan_kernel(u)+epan_kernel(uxl)+epan_kernel(uxu))/h/normal[j];
        }
    }
    // apply contribution to lower ranks
    for (int i = 0; i<global_lower_bound; i++){
        const double x0 = global_x_position[i];
        // fetch local contribution
        for (int j = 0; j<local_size; j++){
            const double h = band_width[j][0];
            const double x = position[j][0];
            const double u = (x0-x)/h;
            const double uxl = ((x-global_min_x)+(x0-global_min_x))/h;
            const double uxu = ((global_max_x-x)+(global_max_x-x0))/h;
            global_contribution[i] += distribution[j]*(epan_kernel(u)+epan_kernel(uxl)+epan_kernel(uxu))/h/normal[j];
        }
    }
    // apply contribution to upper ranks
    for (int i = global_upper_bound; i<size; i++){
        const double x0 = global_x_position[i];
        // fetch local contribution
        for (int j = 0; j<local_size; j++){
            const double h = band_width[j][0];
            const double x = position[j][0];
            const double u = (x0-x)/h;
            const double uxl = ((x-global_min_x)+(x0-global_min_x))/h;
            const double uxu = ((global_max_x-x)+(global_max_x-x0))/h;
            global_contribution[i] += distribution[j]*(epan_kernel(u)+epan_kernel(uxl)+epan_kernel(uxu))/h/normal[j];
        }
    }

    if(domain_decomposed){
        // accumulate global contribution
        rtt_c4::global_sum(global_contribution.data(), size);

        for (int i = 0; i<local_size; i++)
            result[i] += global_contribution[i+global_lower_bound];

    }
    
    return result;
    
}

} // end namespace  kde

#endif // kde_kde_t_hh

//------------------------------------------------------------------------------------------------//
// end of kde/kde.t.hh
//------------------------------------------------------------------------------------------------//

