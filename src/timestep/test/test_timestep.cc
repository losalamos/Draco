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
#include "timestep/test/dummy_package.hh"
#include "timestep/test/test_utils.hh"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

test_timestep::test_timestep()
{
//empty
}

test_timestep::~test_timestep()
{
//empty
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
    dummy_package xxx(mngr);

// Set up a informational advisor to
// contain the current time-step for reference.
// Activating this controller can also be used to
// freeze the time-step at the current value.

    SP<fixed_ts_advisor> sp_dt 
	(new fixed_ts_advisor("Current Time-Step",
			      ts_advisor::req, dt, false));
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

// Now that all the advisors have been set up, perform time cycles

    for (int i=icycle_first; i != icycle_last+1; i++)
    {

	time = time + dt; //end of cycle time
	mngr.set_cycle_data(dt,i,time);

    // Make any user directed changes to controllers

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

    // Pass in the advisors owned by package_XXX for
    // that package to update

   	xxx.advance_state();

    //Compute a new time-step and print results to screen

	dt = mngr.compute_new_timestep();
	mngr.print_summary();
	    
    }

// Dump a list of the advisors to the screen

    mngr.print_advisors();

// Dump the advisor states for visual examination.

    mngr.print_adv_states();

// Confirm that at least some of the output is correct.

    int nd = 5;
    bool passed = mngr.get_cycle() == icycle_last
	&& compare_reals(3.345679e+00, mngr.get_time(),nd)
	&& compare_reals(1.234568e+00, mngr.get_dt(),nd)
	&& compare_reals(1.371742e+00, mngr.get_dt_new(),nd)
	&& mngr.get_controlling_advisor() == "Electron Temperature"
	&& compare_reals(1.000000e-06,sp_min->get_dt_rec(mngr),nd)
	&& compare_reals(9.876543e-01,sp_llr->get_dt_rec(mngr),nd)
	&& compare_reals(1.000000e+00,sp_ovr->get_dt_rec(mngr),nd)
	&& compare_reals(1.234568e+00, sp_dt->get_dt_rec(mngr),nd)
	&& compare_reals(1.481481e+00,sp_ulr->get_dt_rec(mngr),nd)
	&& compare_reals(6.654321e+00, sp_gd->get_dt_rec(mngr),nd)
	&& compare_reals(1.000000e+05,sp_max->get_dt_rec(mngr),nd)
	&& xxx.tests_passed();
	
// Print the status of the test.

    cout << endl;
    cout <<     "*********************************************" << endl;
    if (passed) 
    {
	cout << "**** Time-step Manager Self Test: PASSED ****" << endl;
    }
    else
    {
	cout << "**** Time-step Manager Self Test: FAILED ****" << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;
}


//---------------------------------------------------------------------------//
//                              end of test_timestep.cc
//---------------------------------------------------------------------------//
