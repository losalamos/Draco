//----------------------------------*-C++-*----------------------------------//
// field_ts_advisor.t.cc
// John McGhee
// Mon Aug 24 07:48:00 1998
//---------------------------------------------------------------------------//
// @> Contains the template methods for the field ts_advisor class.
//---------------------------------------------------------------------------//

#include "timestep/field_ts_advisor.hh"

#include "timestep/ts_manager.hh"

#include "ds++/Assert.hh"

#include <algorithm>

#include <stdexcept>

#include <cmath>

using namespace rtt_timestep;

template <class FT>
void field_ts_advisor::set_floor(const FT &y1, double frac)
{
    Require(invariant_satisfied());
    Require(frac > 0.);
    double x1 = -large();
    for (FT::const_iterator py1 = y1.begin(); py1 != y1.end(); py1++) 
    {
	if (*py1 > x1)
	{
	    x1 = *py1;
	}
    }
    x1 = x1*frac;
    if (x1 <= small())
    {
	x1 = small();
    }
    floor_value = x1;
    Ensure(invariant_satisfied());
}


template <class FT>
void field_ts_advisor::update_tstep(const ts_manager &tsm,
				    const FT &q_old, 
				    const FT &q_new) 
{
    Require(invariant_satisfied());
    Require(tsm.get_dt() > 0.0);
//    Require(FT::conformal(q_old,q_new));

    double x1 = 0.;
    double x2 = 0.;
    if (update_method == inf_norm) 
    {
	x2 = 1.;
    }

    FT::const_iterator pq_new = q_new.begin();
    for (FT::const_iterator pq_old = q_old.begin(); 
	 pq_old != q_old.end(); pq_old++,pq_new++) 
    {
	
	if (*pq_old > floor_value && *pq_new > floor_value)
	{
	    double delta_q = std::abs(*pq_new-*pq_old);
	    double q_norm  = *pq_old;
	    double alpha   = delta_q/q_norm;

	    if (alpha < eps()) // Set noise to a hard zero
	    {
		alpha = 0.;
		delta_q = 0.;
	    }
	    if ( update_method == inf_norm )
	    {
		if (alpha > x1) x1=alpha ;
	    }
	    else if ( update_method == a_mean )
	    {
		x1 = x1 + alpha;
		x2 = x2 + 1.;
	    }
	    else if ( update_method == q_mean )
	    {
		x1 = x1 + delta_q;
		x2 = x2 + q_norm;
	    }
	    else if ( update_method == rc_mean )
	    {
		x1 = x1 + alpha*alpha;
		x2 = x2 + alpha;
	    }
	    else if ( update_method == rcq_mean )
	    {
		x1 = x1 + alpha*delta_q;
		x2 = x2 + delta_q;
	    }
	    else
	    {
		throw std::runtime_error("Unrecognized update method flag");
	    }
	}
    }

    if (x1 < small()) 
    {
	dt_rec = large();
    }
    else 
    {
	double fact = x2*fc_value/x1;
	if (fact < small())
	{
	    dt_rec = small();
	}
	else
	{
	    dt_rec = std::max(small(),fact*tsm.get_dt());
	}
    }

    cycle_at_last_update = tsm.get_cycle();
    Ensure(invariant_satisfied());

}


//---------------------------------------------------------------------------//
//                              end of field_ts_advisor.t.cc
//---------------------------------------------------------------------------//


