//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   compton2/Compton2.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Implementation file for compton interface
 * \note   Copyright (C) 2017-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

// headers provided in draco:
#include "compton2/Compton2.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"

#include<iostream>
#include<algorithm>

using UINT = size_t;
using FP = double;
using vec = std::vector<FP>;

namespace rtt_compton2 {

Compton2::Compton2(std::string filename):
    num_temperatures_(0U),
    num_groups_(0U),
    num_leg_moments_(0U),
    num_evals_(0U),
    Ts_(0U),
    Egs_(0U),
    first_groups_(0U),
    indexes_(0U),
    data_(0U),
    derivs_(0U)
{
    Require(filename.length() > 0U);
    std::cout << "Constructor\n";

    // TODO: get rank
    int rank = 0;
    if (rank == 0)
    {
        read_binary(filename);
        broadcast_MPI();
    }
    else
    {
        receive_MPI();
    }

    Ensure(check_class_invariants());
}

void Compton2::broadcast_MPI() const
{
    // TODO: implement
}

void Compton2::receive_MPI()
{
    // TODO: implement
}

void Compton2::read_binary(std::string filename)
{
    // TODO: implement
}

void Compton2::interp_matvec(vec & x,
                             const vec & leftscale,
                             const vec & rightscale,
                             double Te_keV,
                             bool zeroth_moment_only) const
{
    // TODO: Redo interface? Not in-place??

    const vec & L = leftscale;
    const vec & R = rightscale;
    FP const T = std::min(Ts_[num_temperatures_-1], std::max(Ts_[0], Te_keV));
    UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
    // TODO: implement
}

void Compton2::interp_matvec_transpose(vec & xT,
                             const vec & leftscale,
                             const vec & rightscale,
                             double Te_keV,
                             bool zeroth_moment_only) const
{
    // TODO: Redo interface? Not in-place??

    const vec & L = leftscale;
    const vec & R = rightscale;
    FP const T = std::min(Ts_[num_temperatures_-1], std::max(Ts_[0], Te_keV));
    UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
    // TODO: implement
}

void Compton2::interp_sparse_inscat(Sparse_Compton_Matrix & inscat,
                          const vec & leftscale,
                          const vec & rightscale,
                          double Te_keV,
                          bool zeroth_moment_only) const
{
    // TODO: Allow offset to fill in data?

    const vec & L = leftscale;
    const vec & R = rightscale;
    FP const T = std::min(Ts_[num_temperatures_-1], std::max(Ts_[0], Te_keV));
    UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
    // TODO: implement
}

void Compton2::interp_dense_inscat(vec & inscat,
                         const vec & leftscale,
                         const vec & rightscale,
                         double Te_keV,
                         bool zeroth_moment_only) const
{
    const vec & L = leftscale;
    const vec & R = rightscale;
    FP const T = std::min(Ts_[num_temperatures_-1], std::max(Ts_[0], Te_keV));
    UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;

    inscat.resize(num_groups_ * num_groups_ * end_leg) ;
    std::fill(inscat.begin(), inscat.end(), 0.0);

    // TODO: implement
}

void Compton2::interp_linear_outscat(vec & outscat,
                                     const vec & leftscale,
                                     const vec & rightscale,
                                     double Te_keV) const
{
    outscat.resize(num_groups_);
    std::fill(outscat.begin(), outscat.end(), 0.0);

    const vec & L = leftscale;
    const vec & R = rightscale;
    FP const T = std::min(Ts_[num_temperatures_-1], std::max(Ts_[0], Te_keV));
    // TODO: implement
}

void Compton2::interp_nonlinear_diff(vec & nldiff,
                         const vec & leftscale,
                         const vec & rightscale,
                         const vec & flux,
                         double flux_scale,
                         double Te_keV) const
{
    // Need "const" that changes based on DBC level?
    // Require(leftscale.len() == num_groups_);
    // Require(rightscale.len() == num_groups_);
    // Require(flux.len() == num_groups_);
    // Require(flux_scale > 0.0);
    // Require(Te_keV >= 0.0);

    nldiff.resize(num_groups_);
    std::fill(nldiff.begin(), nldiff.end(), 0.0);

    const vec & L = leftscale;
    const vec & R = rightscale;
    FP const fscale = flux_scale;
    FP const T = std::min(Ts_[num_temperatures_-1], std::max(Ts_[0], Te_keV));
    // TODO: implement

    // Ensure(nldiff.len() == num_groups_);
}


} // namespace rtt_compton2

//----------------------------------------------------------------------------//
// End compton2/Compton2.cc
//----------------------------------------------------------------------------//
