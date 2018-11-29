/*----------------------------------*-C++-*----------------------------------//
 * \file   field_ts_advisor_pt.cc
 * \author John McGhee
 * \date   Fri May  1 09:51:28 1998
 * \brief  Explicit template instantiation for the time-step manager test
 *         facility. 
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "field_ts_advisor.t.hh"
#include <vector>

namespace rtt_timestep {

//---------------------------------------------------------------------------//
// Explicit instatiation for FT == vector<double>.
//---------------------------------------------------------------------------//

template void
field_ts_advisor::set_floor<std::vector<double>>(std::vector<double> const &y1,
                                                 double frac);

template void field_ts_advisor::update_tstep<std::vector<double>>(
    ts_manager const &tsm, std::vector<double> const &y1,
    std::vector<double> const &y2);

} // end namespace rtt_timestep

//---------------------------------------------------------------------------//
// end of field_ts_advisor_pt.cc
//---------------------------------------------------------------------------//
