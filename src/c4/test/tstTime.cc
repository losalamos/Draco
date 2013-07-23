//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstTime.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:19:16 2002
 * \brief  Test timing functions in C4.
 * \note   Copyright (C) 2002-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ParallelUnitTest.hh"
#include "../Global_Timer.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>
#include <vector>

// helper macros.
#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define ITFAILS    ut.failure( __LINE__, __FILE__ )

rtt_c4::Global_Timer do_timer("do_global_timer"), do_not_timer("do_not_global_timer");

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void wall_clock_test( rtt_dsxx::UnitTest &ut )
{
    using std::endl;
    using std::cout;
    using std::ostringstream;
    using std::set;
    using std::string;

    using rtt_dsxx::soft_equiv; 
    using rtt_c4::wall_clock_time;
    using rtt_c4::wall_clock_resolution;
    using rtt_c4::Timer;
    using rtt_c4::Global_Timer;

    Global_Timer do_timer("do_timer"), do_not_timer("do_not_timer");
    set<string> active_timers;
    active_timers.insert("do_timer");
    active_timers.insert("do_global_timer");
    Global_Timer::set_global_activity(active_timers);

    double const wcr( rtt_c4::wall_clock_resolution() );
    // double const wcrDeprecated( C4::Wtick() );
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
    // double tolerance(1.0e-14);
    // if( rtt_dsxx::soft_equiv( wcr, wcrDeprecated, tolerance ) )
    // {
    //     ostringstream msg;
    //     msg << "Wtick() and wall_clock_resolution() returned equivalent "
    //         << "values (tolerance = " << tolerance << ")." << endl;
    //     PASSMSG(msg.str());
    // }
    // else
    // {
    //     ostringstream msg;
    //     msg << "The function wall_clock_resolution() returned a value of "
    //         << wcr << " ticks, but the equivalent deprecated function "
    //         << "Wtick() returned a value of " << wcrDeprecated
    //         << " ticks.  These values are not equivalent as they should be."
    //         << endl;
    //     FAILMSG(msg.str());
    // }
    
    Timer t;

    double const prec( 1.8*t.posix_err() );
    
    double begin           = rtt_c4::wall_clock_time();
    //double beginDeprecated = C4::Wtime();

    // if( rtt_dsxx::soft_equiv(begin,beginDeprecated,prec) )
    // {
    //     PASSMSG("Wtime() and wall_clock_time() returned equivalent values.");
    // }
    // else
    // {
    //     ostringstream msg;
    //     msg << "Wtime() and wall_clock_time() failed to return "
    //         << "equivalent values.";
    //     cout.precision(14);
    //     if( beginDeprecated < begin )
    //         msg << "\n\tFound begin < beginDeprecated."
    //             << "\n\tbegin           = " << begin
    //             << "\n\tbeginDeprecated = " << beginDeprecated;
    //     else
    //         msg << "\n\tFound begin != beginDeprecated."
    //             << "\n\tbegin           = " << begin
    //             << "\n\tbeginDeprecated = " << beginDeprecated;
    //     msg << endl;
    //     FAILMSG(msg.str());
    // }
    t.start();

    // do some work
    std::cout << "\nDoing some work..." << std::endl;
    size_t len(20000000);
    std::vector<int> foo(len);
    for( size_t i = 0; i < len; ++i )
        foo[i] = i*3;

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
    
    double const deltaWallTime( t.wall_clock() - 
                                ( t.system_cpu() + t.user_cpu() ) );
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

    cout << "Demonstration of the print() member function via the\n"
         << "\toperator<<(ostream&,Timer&) overloaded operator.\n"
         << endl;

    cout << "Timer = " << t << endl;
        
    //------------------------------------------------------//
    // Do a second timing:
    //------------------------------------------------------//

    cout << "\nCreate a Timer Report after two timing cycles:\n"
         << endl;
    
    t.start();
    for( size_t i = 0; i < len; ++i )
        foo[i]=i*4;
    t.stop();

    t.print( cout, 6 );

    // Test the single line printout
    if( rtt_c4::node() == 0 )
    {
        std::ostringstream timingsingleline;
        t.printline( timingsingleline, 4, 8 );
        // std::cout << "\"" <<  timingsingleline.str() << ".\"" << std::endl;
        // std::cout << "len = " << timingsingleline.str().length() << std::endl;
#ifdef HAVE_PAPI
        if( timingsingleline.str().length() == 42 )
#else
        if( timingsingleline.str().length() == 26 )
#endif
            PASSMSG( "printline() returned a single line of the expected length." );
        else
            FAILMSG( "printline() did not return a line of the expected length." );
    }
    
    //------------------------------------------------------//
    // Check the number of intervals
    //------------------------------------------------------//

    int const expectedNumberOfIntervals(2);
    if( t.intervals() == expectedNumberOfIntervals )
        PASSMSG("Found the expected number of intervals.");
    else
        FAILMSG("Did not find the expected number of intervals.");
          
    //------------------------------------------------------//
    // Check the merge method
    //------------------------------------------------------//

    double old_wall_time = t.sum_wall_clock();
    double old_system_time = t.sum_system_cpu();
    double old_user_time = t.sum_user_cpu();
    double old_intervals = t.intervals();
    t.merge(t);
   
    if (2*old_wall_time==t.sum_wall_clock() &&
        2*old_system_time==t.sum_system_cpu() &&
        2*old_user_time==t.sum_user_cpu() &&
        2*old_intervals==t.intervals())
        PASSMSG("merge okay");
    else
        FAILMSG("merge NOT okay");

    //------------------------------------------------------------//
    // Check PAPI data
    //------------------------------------------------------------//

    long long cachemisses = t.sum_L2_cache_misses();
    long long cachehits   = t.sum_L2_cache_hits();
    long long flops       = t.sum_floating_operations();

#ifdef HAVE_PAPI
    
    std::cout << "PAPI metrics report:\n"
              << "   Cache misses : " << cachemisses << "\n"
              << "   Cache hits   : " << cachehits << "\n"
              << "   FLOP         : " << flops << std::endl;

    if( cachemisses == 0 && cachehits == 0  && flops == 0 )
        FAILMSG( "PAPI metrics returned 0 when PAPI was available.");
    else
        PASSMSG( "PAPI metrics returned >0 values when PAPI is available.");
    
#else
    if( cachemisses == 0 && cachehits == 0  && flops == 0 )
        PASSMSG( "PAPI metrics return 0 when PAPI is not available.");
    else
        FAILMSG( "PAPI metrics did not return 0 when PAPI was not available.");
#endif
    
    return;
}

//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {   // UNIT TESTS
      	wall_clock_test(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstTime, " 
                  << err.what() << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstTime, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstTime.cc
//---------------------------------------------------------------------------//
