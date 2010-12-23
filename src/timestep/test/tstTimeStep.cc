//----------------------------------*-C++-*----------------------------------//
/*! \file   tstTimeStep.c
 *  \author John McGhee
 *  \date   Fri May  1 09:43:49 1998
 *  \brief  A driver for the time-step manager test facility.
 *  \note   Copyright (C) 1998-2010 Los Alamos National Security, LLC.
 *          All rights reserved.  */
//---------------------------------------------------------------------------//
//! \version $Id$
//---------------------------------------------------------------------------//

#include "dummy_package.hh"
#include "timestep_test.hh"
#include "../ts_manager.hh"
#include "../Release.hh"
#include "../fixed_ts_advisor.hh"
#include "../ratio_ts_advisor.hh"
#include "../target_ts_advisor.hh"
#include "../field_ts_advisor.hh"

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

// forward declaration
void run_tests();
void check_field_ts_advisor();

// Main program
int main ( int argc, char *argv[] )
{
    using namespace std;
    using namespace rtt_c4;

    initialize(argc, argv);

    // version tag
    if( node() == 0 )
        cout << argv[0] << ": version " << rtt_timestep::release() << endl;
    for( int arg=1; arg < argc; arg++ )
	if( string(argv[arg]) == "--version" )
	{ finalize(); return 0; }
	
    try
    { 	// Run the tests...
	run_tests();
        check_field_ts_advisor();
    }
    catch( const std::exception &err )
    {
	cerr << "exception: " << err.what() << endl;
        rtt_c4::abort();
	return 1;
    }
    catch( ... )
    {
	std::cerr << "unknown exception" << std::endl;
        rtt_c4::abort();
	return 1;
    }

    // Report result
    {
	HTSyncSpinLock slock;
	
	// status of test
	cout <<     "\n*********************************************";
	if( rtt_timestep_test::passed ) 
	    cout << "\n**** timestep Test: PASSED on " << node();
	else
	    cout << "\n**** timestep Test: FAILED on " << node();
	cout <<     "\n*********************************************\n\n";
    }
    
    global_barrier();    
    cout << "Done testing tstTime on " << node() << endl;
    finalize();
    return 0;
}

//---------------------------------------------------------------------------//
// Actual tests go here.
//---------------------------------------------------------------------------//

void run_tests()
{
    using std::cout;
    using std::endl;
    using rtt_dsxx::SP;
    using rtt_dsxx::soft_equiv;

    using rtt_timestep::fixed_ts_advisor;
    using rtt_timestep::ratio_ts_advisor;
    using rtt_timestep::target_ts_advisor;
    using rtt_timestep::ts_manager;
    using rtt_timestep::ts_advisor;
    
    using rtt_timestep_test::dummy_package;

    // Initial values;
    double graphics_time( 10.0 );
    double dt_min(         0.000001 );
    double dt_max(    100000.0 );
    double override_dt(    1.0 );
    double dt(             1.0 );
    double time(           0.  );

    int icycle_first( 1 );
    int icycle_last(  3 );

    bool override_flag( false );

    ts_manager mngr;
    dummy_package xxx( mngr );

    // Set up a informational advisor to contain the current time-step for
    // reference.  Activating this controller can also be used to freeze the
    // time-step at the current value.

    SP< fixed_ts_advisor > sp_dt(
	new fixed_ts_advisor( "Current Time-Step",
			      ts_advisor::req, dt, false) );
    mngr.add_advisor( sp_dt );

    // Set up a required time-step to be activated at the user's discretion
    
    SP< fixed_ts_advisor > sp_ovr(
	new fixed_ts_advisor( "User Override",
			      ts_advisor::req, 
			      override_dt, false) );
    mngr.add_advisor( sp_ovr );
    
    // Set up a min timestep

    SP< fixed_ts_advisor > sp_min(
	new fixed_ts_advisor( "Minimum",
			      ts_advisor::min, 
			      ts_advisor::ts_small()) );
    mngr.add_advisor( sp_min );
    sp_min -> set_fixed_value( dt_min );

    // Set up a lower limit on the timestep rate of change

    SP< ratio_ts_advisor > sp_llr(
	new ratio_ts_advisor( "Rate of Change Lower Limit",
			       ts_advisor::min, 0.8 ) );
    mngr.add_advisor( sp_llr );

    // Test the accessor
    {
        double tmp( sp_llr->get_ratio() );
        Check( rtt_dsxx::soft_equiv( tmp, 0.8 ) );
    }

    // Set up an upper limit on the time-step rate of change

    SP< ratio_ts_advisor > sp_ulr(
	new ratio_ts_advisor("Rate of Change Upper Limit") );
    mngr.add_advisor( sp_ulr );

    // Set up an advisor to watch for an upper limit on the time-step.
    
    SP< fixed_ts_advisor > sp_max( new fixed_ts_advisor( "Maximum" ) );
    mngr.add_advisor( sp_max );
    sp_max -> set_fixed_value( dt_max );
    
    // Set up a target time advisor

    SP< target_ts_advisor > sp_gd(
	new target_ts_advisor( "Graphics Dump",
			       ts_advisor::max, graphics_time) );
    mngr.add_advisor( sp_gd );

    // Now that all the advisors have been set up, perform time cycles

    for( int i=icycle_first; i <= icycle_last; ++i )
    {

	time = time + dt; //end of cycle time
	mngr.set_cycle_data( dt, i, time );

	// Make any user directed changes to controllers

	sp_dt -> set_fixed_value( dt );
	if( override_flag )
	{
	    sp_ovr -> activate();
	    sp_ovr -> set_fixed_value( override_dt );
	}
	else
	{
	    sp_ovr -> deactivate();
	}
	
	// Pass in the advisors owned by package_XXX for
	// that package to update

   	xxx.advance_state();
	
	// Compute a new time-step and print results to screen
	
	dt = mngr.compute_new_timestep();
	mngr.print_summary();
    }
    
    // Dump a list of the advisors to the screen
    
    mngr.print_advisors();
    
    // Dump the advisor states for visual examination.
    
    mngr.print_adv_states();
    
    //------------------------------------------------------------//
    // Confirm that at least some of the output is correct.
    //------------------------------------------------------------//
    
    // Reference Values:
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
    double const ref11( 1.371742 );
    double const ref12( 1.496914 );
    double const ref13( 2.716049 );

    // Check final values:
    // ------------------------------
    if( mngr.get_cycle() == icycle_last )
    {
	PASSMSG("get_cycle() returned the expected cycle index.");
    }
    else
    {
	FAILMSG("get_cycle() failed to return the expected cycle index.");
    }

    if( mngr.get_controlling_advisor() == "Electron Temperature" ) 
    {
	PASSMSG("get_controlling_advisor() returned expected string.");
    }
    else
    {
	FAILMSG("get_controlling_advisor() failed to return the expected string.");
    }

    if( soft_equiv( ref11, xxx.get_dt_rec_te(), prec ) )
    {
	PASSMSG("get_dt_rec_te() gave expected value.");
    }
    else
    {
	FAILMSG("get_dt_rec_te() did not give expected value.");
    }

    if( soft_equiv( ref12, xxx.get_dt_rec_ti(), prec ) )
    {
	PASSMSG("get_dt_rec_ti() gave expected value.");
    }
    else
    {
	FAILMSG("get_dt_rec_ti() did not give expected value.");
    }

    if( soft_equiv( ref13, xxx.get_dt_rec_ri(), prec ) )
    {
	PASSMSG("get_dt_rec_ri() gave expected value.");
    }
    else
    {
	FAILMSG("get_dt_rec_ri() did not give expected value.");
    }

    if( soft_equiv( ref1, mngr.get_time(), prec ) )
    {
	PASSMSG("get_time() gave expected value.");
    }
    else
    {
	FAILMSG("get_time() did not give expected value.");
    }

    if( soft_equiv( ref2, mngr.get_dt(), prec ) )
    {
	PASSMSG("get_dt() gave expected value.");
    }
    else
    {
	FAILMSG("get_dt() did not give expected value.");
    }

    if( soft_equiv( ref3, mngr.get_dt_new(), prec ) )
    {
	PASSMSG("get_dt_new() gave expected value.");
    }
    else
    {
	FAILMSG("get_dt_new() did not give expected value.");
    }

    if( soft_equiv( ref4, sp_min->get_dt_rec(mngr), prec ) )
    {
	PASSMSG(" sp_min->get_dt_rec(mngr) gave expected value.");
    }
    else
    {
	FAILMSG( "sp_min->get_dt_rec(mngr) did not give expected value.");
    }

    if( soft_equiv( ref5, sp_llr->get_dt_rec(mngr), prec ) )
    {
	PASSMSG(" sp_llr->get_dt_rec(mngr) gave expected value.");
    }
    else
    {
	FAILMSG( "sp_llr->get_dt_rec(mngr) did not give expected value.");
    }

    if( soft_equiv( ref6, sp_ovr->get_dt_rec(mngr), prec ) )
    {
	PASSMSG(" sp_ovr->get_dt_rec(mngr) gave expected value.");
    }
    else
    {
	FAILMSG( "sp_ovr->get_dt_rec(mngr) did not give expected value.");
    }

    if( soft_equiv( ref7, sp_dt->get_dt_rec(mngr), prec ) )
    {
	PASSMSG(" sp_dt->get_dt_rec(mngr) gave expected value.");
    }
    else
    {
	FAILMSG( "sp_dt->get_dt_rec(mngr) did not give expected value.");
    }

    if( soft_equiv( ref8, sp_ulr->get_dt_rec(mngr), prec ) )
    {
	PASSMSG(" sp_ulr->get_dt_rec(mngr) gave expected value.");
    }
    else
    {
	FAILMSG( "sp_ulr->get_dt_rec(mngr) did not give expected value.");
    }

    if( soft_equiv( ref9, sp_gd->get_dt_rec(mngr), prec ) )
    {
	PASSMSG(" sp_gd->get_dt_rec(mngr) gave expected value.");
    }
    else
    {
	FAILMSG( "sp_gd->get_dt_rec(mngr) did not give expected value.");
    }

    if( soft_equiv( ref10, sp_max->get_dt_rec(mngr), prec ) )
    {
	PASSMSG(" sp_max->get_dt_rec(mngr) gave expected value.");
    }
    else
    {
	FAILMSG( "sp_max->get_dt_rec(mngr) did not give expected value.");
    }

    // Check to make sure all processes passed.
    
    int npassed = rtt_timestep_test::passed ? 1 : 0;
    C4::gsum( npassed );

    if( npassed == C4::nodes() )
    {
	PASSMSG("All tests passed on all procs.");
    }
    else
    {
	std::ostringstream msg;
	msg << "Some tests failed on processor " << C4::node() << std::endl;
	FAILMSG( msg.str() );
    }

    return;
}

void check_field_ts_advisor()
{
    std::cout << "\nChecking the field_ts_advisor class...\n" << std::endl;
    
    rtt_timestep::field_ts_advisor ftsa;

    // Check manipulators
    std::cout << "Setting Frac Change to 1.0..." << std::endl;
    ftsa.set_fc( 1.0 );
    std::cout << "Setting Floor Value to 0.0001..." << std::endl;
    ftsa.set_floor( 0.0001 );
    std::cout << "Setting Update Method to q_mean..." << std::endl;
    ftsa.set_update_method( rtt_timestep::field_ts_advisor::q_mean );

    // Dump the state to an internal buffer and inspect the results
    std::ostringstream msg;
    ftsa.print_state( msg );
    std::cout << msg.str() << std::endl;

    { // Check the Fraction Change value
        std::string const expected( "Fract Change   : 1" );

        // find the line of interest
        std::string output( msg.str() );
        size_t beg( output.find("Fract Change") );
        if( beg == std::string::npos )
        {
            FAILMSG("Did not find expected string!");
            return;
        }
        size_t end( output.find_first_of("\n",beg) );
        if( beg == std::string::npos )
        {
            FAILMSG("Did not find expected string!");
            return;
        }
        std::string line( output.substr(beg,end-beg) );
        if( line == expected ) {
            PASSMSG("'Fract Change' was set correctly."); }
        else {
            FAILMSG("Failed to set 'Fract Change' correctly."); }
    }
    
   { // Check the Floor value
        std::string const expected( "Floor Value    : 0.0001" );

        // find the line of interest
        std::string output( msg.str() );
        size_t beg( output.find("Floor Value") );
        if( beg == std::string::npos )
        {
            FAILMSG("Did not find expected string!");
            return;
        }
        size_t end( output.find_first_of("\n",beg) );
        if( beg == std::string::npos )
        {
            FAILMSG("Did not find expected string!");
            return;
        }
        std::string line( output.substr(beg,end-beg) );
        if( line == expected ) {
            PASSMSG("'Floor Value' was set correctly."); }
        else {
            FAILMSG("Failed to set 'Floor Value' correctly."); }
    }

   { // Check the Update Method value
        std::string const expected( "Update Method  : weighted by field value" );

        // find the line of interest
        std::string output( msg.str() );
        size_t beg( output.find("Update Method") );
        if( beg == std::string::npos )
        {
            FAILMSG("Did not find expected string!");
            return;
        }
        size_t end( output.find_first_of("\n",beg) );
        if( beg == std::string::npos )
        {
            FAILMSG("Did not find expected string!");
            return;
        }
        std::string line( output.substr(beg,end-beg) );
        if( line == expected ) {
            PASSMSG("'Update Method' was set correctly."); }
        else {
            FAILMSG("Failed to set 'Update Method' correctly."); }
    }
      
    return;    
}


//---------------------------------------------------------------------------//
//                         end of tstTimeStep.c
//---------------------------------------------------------------------------//
