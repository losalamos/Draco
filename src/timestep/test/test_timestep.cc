//----------------------------------*-C++-*----------------------------------//
// test_timestep.cc
// John McGhee
// Fri May  1 09:44:46 1998
//---------------------------------------------------------------------------//
// @> Performs test of time-step manager facility.
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <cmath>
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "../ts_manager.hh"
#include "../fixed_ts_advisor.hh"
#include "../ratio_ts_advisor.hh"
#include "../target_ts_advisor.hh"
#include "../test/dummy_package.hh"
#include "test_timestep.hh"

using namespace rtt_timestep;

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
    using std::cout;
    using std::endl;
    using rtt_dsxx::SP;
    using rtt_dsxx::soft_equiv;

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

    double const prec( 1.0e-5 );
    double const ref1( 3.345679 );
    double const ref2( 1.234568 );
    double const ref3( 1.371742 );
    double const ref4( 1.000000e-06 );
    double const ref5( 9.876543e-01 );
    double const ref6( 1.0 );
    double const ref7( 1.234568 );
    double const ref8( 1.481481 );
    double const ref9( 6.654321 );
    double const ref10( 1.000000e+05 );

    bool passed = mngr.get_cycle() == icycle_last
 	&& mngr.get_controlling_advisor() == "Electron Temperature"
 	&& xxx.tests_passed()
	&& soft_equiv( ref1, mngr.get_time(), prec )
	&& soft_equiv( ref2, mngr.get_dt(), prec )
	&& soft_equiv( ref3, mngr.get_dt_new(), prec )
	&& soft_equiv( ref4, sp_min->get_dt_rec(mngr), prec )
	&& soft_equiv( ref5, sp_llr->get_dt_rec(mngr), prec )
	&& soft_equiv( ref6, sp_ovr->get_dt_rec(mngr), prec )
 	&& soft_equiv( ref7, sp_dt->get_dt_rec(mngr),prec)
 	&& soft_equiv( ref8, sp_ulr->get_dt_rec(mngr),prec)
 	&& soft_equiv( ref9, sp_gd->get_dt_rec(mngr),prec)
 	&& soft_equiv( ref10,sp_max->get_dt_rec(mngr),prec);

    // Check to make sure all processes passed.
    
    int npassed = passed ? 1 : 0;
    C4::gsum(npassed);

    passed = npassed == C4::nodes();
    
// Print the status of the test.

    if (C4::node() == 0)
    {
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
}


//---------------------------------------------------------------------------//
//                              end of test_timestep.cc
//---------------------------------------------------------------------------//
