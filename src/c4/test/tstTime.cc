//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstTime.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:19:16 2002
 * \brief  Test timing functions in C4.
 * \note   Copyright (C) 2002-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Release.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "../Timer.hh"
#include "c4_test.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void wall_clock_test()
{
    using std::endl;
    using std::cout;
    using std::ostringstream;

    using rtt_dsxx::soft_equiv; 
    using rtt_c4::wall_clock_time;
    using rtt_c4::wall_clock_resolution;
    using rtt_c4::Timer;

    double const wcr( rtt_c4::wall_clock_resolution() );
    double const wcrDeprecated( C4::Wtick() );
    if( wcr > 0.0 && wcr <= 100.0)
    {
        ostringstream msg;
        msg << "The timer has a wall clock resoution of "
            << wcr << " ticks." << endl;
        PASSMSG(msg.str());
    }
    else
    {
        ostringstream msg;
        msg << "The timer does not appear to have a reasonable resolution."
            << " rtt_c4::wall_clock_resolution() = " << wcr << " ticks."
            << endl;
        FAILMSG(msg.str());
    }
    double tolerance(1.0e-14);
    if( rtt_dsxx::soft_equiv( wcr, wcrDeprecated, tolerance ) )
    {
        ostringstream msg;
        msg << "Wtick() and wall_clock_resolution() returned equivalent "
            << "values (tolerance = " << tolerance << ")." << endl;
        PASSMSG(msg.str());
    }
    else
    {
        ostringstream msg;
        msg << "The function wall_clock_resolution() returned a value of "
            << wcr << " ticks, but the equivalent deprecated function "
            << "Wtick() returned a value of " << wcrDeprecated
            << " ticks.  These values are not equivalent as they should be."
            << endl;
        FAILMSG(msg.str());
    }
    
    Timer t;

    double const prec( 1.75*t.posix_err() );
    
    double begin           = rtt_c4::wall_clock_time();
    double beginDeprecated = C4::Wtime();

    if( rtt_dsxx::soft_equiv(begin,beginDeprecated,prec) )
    {
        PASSMSG("Wtime() and wall_clock_time() returned equivalent values.");
    }
    else
    {
        ostringstream msg;
        msg << "Wtime() and wall_clock_time() failed to return "
            << "equivalent values.";
        cout.precision(14);
        if( beginDeprecated < begin )
            msg << "\n\tFound begin < beginDeprecated."
                << "\n\tbegin           = " << begin
                << "\n\tbeginDeprecated = " << beginDeprecated;
        else
            msg << "\n\tFound begin != beginDeprecated."
                << "\n\tbegin           = " << begin
                << "\n\tbeginDeprecated = " << beginDeprecated;
        msg << endl;
        FAILMSG(msg.str());
    }
    t.start();
    
    for( int i = 0; i < 200000000; i++ )
    { /* empty */
    }

    double end = rtt_c4::wall_clock_time();
    t.stop();

    double const error( t.wall_clock() - (end-begin) );
    if( std::fabs(error) <= prec )
    {
	PASSMSG("wall_clock() value looks ok.");
    }
    else
    {
	ostringstream msg;
	msg << "t.wall_clock() value does not match the expected value."
	    << "\n\tend            = " << end
	    << "\n\tbegin          = " << begin
	    << "\n\tend-begin      = " << end - begin
	    << "\n\tt.wall_clock() = " << t.wall_clock()
	    << "\n\tprec           = " << prec << endl;
	FAILMSG(msg.str());
    }

    //---------------------------------------------------------------------//
    // Ensure that system + user <= wall
    //
    // Due to round off errors, the wall clock time might be less than the
    // system + user time.  But this difference should never exceed
    // t.posix_err(). 
    //---------------------------------------------------------------------//
    
    double const deltaWallTime( t.wall_clock() - (
				    t.system_cpu() + t.user_cpu() ) );
#ifdef _MSC_VER
    double const time_resolution( 1.0 );  
#else
    double const time_resolution( prec );  
#endif
    if( deltaWallTime > 0.0 || std::fabs(deltaWallTime) <= time_resolution )
    {
	ostringstream msg;
	msg << "The sum of cpu and user time is less than or equal to the\n\t"
	    << "reported wall clock time (within error bars = " << time_resolution
	    << " secs.)." << endl;
	PASSMSG(msg.str());
    }
    else
    {
	ostringstream msg;
	msg << "The sum of cpu and user time exceeds the reported wall "
	    << "clock time.  Here are the details:"
	    << "\n\tposix_error() = " << prec << " sec."
	    << "\n\tdeltaWallTime = " << deltaWallTime  << " sec."
 	    << "\n\tSystem time   = " << t.system_cpu() << " sec."
 	    << "\n\tUser time     = " << t.user_cpu()   << " sec."
 	    << "\n\tWall time     = " << t.wall_clock() << " sec."
	    << endl;
	FAILMSG(msg.str());
    }

    //------------------------------------------------------//
    // Demonstrate print functions:
    //------------------------------------------------------//

    cout << "\nDemonstration of the print() member function via the\n"
         << "\toperator<<(ostream&,Timer&) overloaded operator.\n"
         << endl;

    cout << "Timer = " << t << endl;
        
    //------------------------------------------------------//
    // Do a second timing:
    //------------------------------------------------------//

    cout << "\nCreate a Timer Report after two timing cycles:\n"
         << endl;
    
    t.start();
    for( int i = 0; i < 200000000; i++ )
    { /* empty */
    }
    t.stop();

    t.print( cout, 6 );

    //------------------------------------------------------//
    // Check the number of intervals
    //------------------------------------------------------//

    int const expectedNumberOfIntervals(2);
    if( t.intervals() == expectedNumberOfIntervals )
        PASSMSG("Found the expected number of intervals.")
    else
        FAILMSG("Did not find the expected number of intervals.")
    
    return;
}

//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    using std::cout;
    using std::endl;
    using std::string;
    
    rtt_c4::initialize( argc, argv );

    // version tag
    if( rtt_c4::node() == 0 )
        cout << argv[0] << ": version " << rtt_c4::release() 
             << endl;

    for( int arg = 1; arg < argc; arg++ )
	if( string( argv[arg] ) == "--version" )
	{
	    rtt_c4::finalize();
	    return 0;
	}

//---------------------------------------------------------------------------//
// UNIT TESTS
//---------------------------------------------------------------------------//
    try 
    { 		
	wall_clock_test();
    }
    catch( rtt_dsxx::assertion &assert )
    {
	cout << "While testing tstTime, " << assert.what()
	     << endl;
	rtt_c4::abort();
	return 1;
    }

//---------------------------------------------------------------------------//
// Print status of test
//---------------------------------------------------------------------------//
    {
	// status of test
	cout <<   "\n*********************************************\n";
	if( rtt_c4_test::passed )
	    cout << "**** tstTime Test: PASSED on " << rtt_c4::node() << endl;
        else
            cout << "**** tstTime Test: FAILED on " << rtt_c4::node() << endl;
	cout <<     "*********************************************\n" << endl;
    }

    rtt_c4::global_barrier();
    cout << "Done testing tstTime on " << rtt_c4::node() << endl;
    rtt_c4::finalize();
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstTime.cc
//---------------------------------------------------------------------------//
