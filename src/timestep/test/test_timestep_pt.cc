//----------------------------------*-C++-*----------------------------------//
// test_timestep_pt.cc
// John McGhee
// Fri May  1 09:51:28 1998
//---------------------------------------------------------------------------//
// @> Template instantiation for the time-step manager test facility.
//---------------------------------------------------------------------------//

#include "../field_ts_advisor.t.hh"
#include <vector>

using namespace rtt_timestep;

typedef std::vector<double> FTVD;

template void field_ts_advisor::set_floor(const FTVD &y1, double frac); 

template void field_ts_advisor::update_tstep(const ts_manager &tsm,
					     const FTVD &y1, const FTVD &y2);

//---------------------------------------------------------------------------//
//                              end of test_timestep_pt.cc
//---------------------------------------------------------------------------//
