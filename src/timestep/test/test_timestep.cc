//----------------------------------*-C++-*----------------------------------//
// test_timestep.cc
// John McGhee
// Fri May  1 09:44:46 1998
//---------------------------------------------------------------------------//
// @> Performs test of time-step manager facility.
//---------------------------------------------------------------------------//

#include "timestep/test/test_timestep.hh"

#include "timestep/ts_manager.hh"
#include "timestep/fixed_ts_advisor.hh"
#include "timestep/ratio_ts_advisor.hh"
#include "timestep/target_ts_advisor.hh"
#include "timestep/field_ts_advisor.hh"

#include <vector>
#include <iostream>
using std::cerr;
using std::endl;

using std::vector;

test_timestep::test_timestep()
{
}

test_timestep::~test_timestep()
{
}

vector<double> operator*(double lhs, const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs*rhs[i];

    return results;
}

void package_XXX(field_ts_advisor &telect_adv, 
		 field_ts_advisor &tion_adv,
		 field_ts_advisor &radi_adv, 
		 const double dt, const int icycle)
{

// Create a set of dummy arrays to serve as control fields for
// use in exercizing the various advisors.

    const double a1[] = {1., 10., 11., 3., 2., 5., 5., 6.7};
    int sizea = sizeof(a1)/sizeof(a1[0]);
    vector<double> te_old(a1,a1+sizea);
    vector<double> te_new = 1.09*te_old;
    vector<double> ti_old=0.97*te_old;
    vector<double> ti_new=1.05*te_old;
    vector<double> ri_old=1.10*te_old;
    vector<double> ri_new=1.15*te_old;

// Set a floor for the electron temperature controller, to
// execcize this method. Just accelpt the default floor on the
// other controllers.

    telect_adv.set_floor(te_new,0.001);

// Get a new time-step from each of the advisors that
// belong to this package.

    telect_adv.update_tstep(te_old, te_new, dt, icycle);
    tion_adv.update_tstep(ti_old, ti_new, dt, icycle);
    radi_adv.update_tstep(ri_old, ri_new, dt, icycle);
}


void test_timestep::execute_test()
{

    using dsxx::SP;

    double graphics_time = 10.;
    double dt_min = 0.000001;
    double dt_max = 100000.;
    bool override_flag = false;
    double override_dt = 1.;
    int icycle_first = 1;
    int icycle_last = 3;
    double dt = 1.;
    double time = 0.;
    ts_manager mngr;

// Set up a informational advisor to
// contain the current time-step for reference

    SP<fixed_ts_advisor> sp_dt 
	(new fixed_ts_advisor("Current Time-Step",
			      ts_advisor::max, dt, false));
    mngr.add_advisor(sp_dt);

// Set up a required time-step to be activated
// at the user's discretion

    SP<fixed_ts_advisor> sp_ovr 
	(new fixed_ts_advisor("User Override",
			      ts_advisor::req, 
			      override_dt, false));
    mngr.add_advisor(sp_ovr);

// Set up a min timestep

    SP<fixed_ts_advisor> sp_min 
	(new fixed_ts_advisor("Minimum",
			      ts_advisor::min, 
			      ts_advisor::small()));
    mngr.add_advisor(sp_min);
    sp_min -> set_fixed_value(dt_min);

// Set up a lower limit on the timestep rate of change

    SP<ratio_ts_advisor> sp_llr 
	(new ratio_ts_advisor( "Rate of Change Lower Limit",
			       ts_advisor::min, 0.8 ) );
    mngr.add_advisor(sp_llr);

// Set up a Electron Temperature advisor

    SP<field_ts_advisor> sp_te 
	(new field_ts_advisor( "Electron Temperature",
			       ts_advisor::max,
			       field_ts_advisor::a_mean));
    mngr.add_advisor(sp_te);

// Set up a Ion-Temperature advisor

    SP<field_ts_advisor> sp_ti 
	(new field_ts_advisor("Ion Temperature",
			      ts_advisor::max,
			      field_ts_advisor::rc_mean));
    mngr.add_advisor(sp_ti);

// Set up a Radiation-Intensity advisor

    SP<field_ts_advisor> sp_ri 
	(new field_ts_advisor("Radiation Intensity",
			      ts_advisor::max,
			      field_ts_advisor::rcq_mean));
    mngr.add_advisor(sp_ri);

// Set up an upper limit on the time-step rate of change

    SP<ratio_ts_advisor> sp_ulr 
	(new ratio_ts_advisor("Rate of Change Upper Limit"));
    mngr.add_advisor(sp_ulr);

// Set up an advisor to watch for an upper limit on the time-step.

    SP<fixed_ts_advisor> sp_max (new fixed_ts_advisor("Maximum"));
    mngr.add_advisor(sp_max);
    sp_max -> set_fixed_value(dt_max);

// Set up a target time advisor

    SP<target_ts_advisor> sp_gd 
	(new target_ts_advisor("Graphics Dump",
			       ts_advisor::max, graphics_time));
    mngr.add_advisor(sp_gd);

// Dump a summary to the screen

    mngr.print_advisors();

// Now that all the advisors have been set up, perform time cycles

    for (int i=icycle_first; i != icycle_last+1; i++)
    {

	time = time + dt; //end of cycle time

    // Update the advisors owned by the host
	sp_dt -> set_fixed_value(dt);
	if (override_flag)
	{
	    sp_ovr -> activate();
	    sp_ovr -> set_fixed_value(override_dt);
	}
	else
	{
	    sp_ovr -> deactivate();
	}
    	sp_dt  -> update_tstep(i);
	sp_ovr -> update_tstep(i);
	sp_min -> update_tstep(i);
	sp_llr -> update_tstep(dt,i);
	sp_ulr -> update_tstep(dt,i);
	sp_max -> update_tstep(i);
	sp_gd  -> update_tstep(time,i);

    // Pass in the advisors owned by package_XXX for
    // that package to update

   	package_XXX(*sp_te, *sp_ti, *sp_ri, dt, i);

    //Compute a new time-step and print results to screen

	dt = mngr.compute_new_timestep(dt,i,time);
	mngr.print_summary();
	    
    }

//Dump some advisor states for examination

    sp_te  -> print_state();
    sp_ulr -> print_state();
    sp_gd  -> print_state();
    sp_max -> print_state();
    sp_min -> print_state();
    sp_dt  -> print_state();
}


//---------------------------------------------------------------------------//
//                              end of test_timestep.cc
//---------------------------------------------------------------------------//
