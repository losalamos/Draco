//----------------------------------*-C++-*----------------------------------//
// test_timestep_pt.cc
// John McGhee
// Fri May  1 09:51:28 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "timestep/field_ts_advisor.t.cc"
#include <vector>

typedef std::vector<double> FTVD;

template void field_ts_advisor::set_floor(const FTVD &y1, double frac); 

template void field_ts_advisor::update_tstep(const FTVD &y1, const FTVD &y2,
					     double current_dt, int cycle_);

//---------------------------------------------------------------------------//
//                              end of test_timestep_pt.cc
//---------------------------------------------------------------------------//
