//----------------------------------*-C++-*----------------------------------//
// tstTimeStep.c
// John McGhee
// Fri May  1 09:43:49 1998
//---------------------------------------------------------------------------//
// @> A driver for the time-step manager test facility.
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"

#include "../ts_manager.hh"
#include "../Release.hh"
#include "../fixed_ts_advisor.hh"
#include "../ratio_ts_advisor.hh"
#include "../target_ts_advisor.hh"

#include "dummy_package.hh"
#include "timestep_test.hh"

// forward declaration
void run_tests();

// Main program
int main ( int argc, char *argv[] )
{
    using std::cout;
    using std::endl;

    C4::Init(argc, argv);

    // version tag
    for( int arg=1; arg < argc; arg++ )
    {
	if( std::string(argv[arg]) == "--version" )
	{
	    if( C4::node() == 0 )
		cout << argv[0] << ": version " 
		     << rtt_timestep::release() << endl;
	    C4::Finalize();
	    return 0;
	}
    }
	
    try
    { 	// Run the tests...
	run_tests();
    }
    catch( const rtt_dsxx::assertion &ass )
    {
	std::cerr << "assert failed: " << ass.what() << std::endl;
	C4::Finalize();
	return 1;
    }
    catch( const std::exception &ass )
    {
	std::cerr << "exception: " << ass.what() << std::endl;
	C4::Finalize();
	return 1;
    }
    catch( ... )
    {
	std::cerr << "unknown exception" << std::endl;
	C4::Finalize();
	return 1;
    }

    // Report result
    {
	C4::HTSyncSpinLock slock;
	
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if( rtt_timestep_test::passed ) 
	{
	    cout << "**** timestep Test: PASSED on " << C4::node() << endl;
	}
	else
	{
	    cout << "**** timestep Test: FAILED on " << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();
    
    cout << "Done testing tstTime on " << C4::node() << endl;
    
    C4::Finalize();
    
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
			      ts_advisor::small()) );
    mngr.add_advisor( sp_min );
    sp_min -> set_fixed_value( dt_min );

    // Set up a lower limit on the timestep rate of change

    SP< ratio_ts_advisor > sp_llr(
	new ratio_ts_advisor( "Rate of Change Lower Limit",
			       ts_advisor::min, 0.8 ) );
    mngr.add_advisor( sp_llr );

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


//---------------------------------------------------------------------------//
//                         end of tstTimeStep.c
//---------------------------------------------------------------------------//
