//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstTime.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:19:16 2002
 * \brief  Test timing functions in C4.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "ds++/Soft_Equivalence.hh"
#include "../Release.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "../Timer.hh"
#include "c4_test.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void wall_clock_test()
{
    using rtt_dsxx::soft_equiv;
    using rtt_c4::wall_clock_time;
    using rtt_c4::wall_clock_resolution;
    using rtt_c4::Timer;
    
    Timer t;

    double begin = rtt_c4::wall_clock_time();
    t.start();
    
    for( int i = 0; i < 200000000; i++ )
    { /* empty */
    }

    double end = rtt_c4::wall_clock_time();
    t.stop();

    double const prec( t.posix_err() );
    if( soft_equiv( end-begin, t.wall_clock(), prec ) )
    {
	PASSMSG("wall_clock() value looks ok.");
    }
    else
    {
	FAILMSG("wall_clock() value does not match expected value.");
    }

    //---------------------------------------------------------------------//
    // Ensure that system + user <= wall
    //
    // Due to round off errors, the wall clock time might be less than the
    // system + user time.  But this difference should never exceed
    // t.posix_err(). 
    //---------------------------------------------------------------------//
    
    double const deltaWallTime( t.wall_clock() - ( t.system_cpu() + t.user_cpu() ) );
    
    if( deltaWallTime > 0.0 || std::fabs(deltaWallTime) < prec )
    {
	std::ostringstream msg;
	msg << "The sum of cpu and user time is less than or equal to the\n\t"
	    << "reported wall clock time (within error bars = " << prec
	    << " secs.)." << std::endl;
	PASSMSG(msg.str());
    }
    else
    {
	std::ostringstream msg;
	msg << "The sum of cpu and user time exceeds the reported wall "
	    << "clock time." << std::endl
	    << "\tSystem time: " << t.system_cpu() << std::endl
	    << "\tUser time  : " << t.user_cpu()   << std::endl
	    << "\tWall time  : " << t.wall_clock() << std::endl;
	FAILMSG(msg.str());
    }

    //------------------------------------------------------//
    // Demonstrate print functions:
    //------------------------------------------------------//

    t.print( std::cout, 6 );

        
    //------------------------------------------------------//
    // Do a second timing:
    //------------------------------------------------------//

    std::cout << "\nCreate a Timer Report after two timing cycles:\n"
	      << std::endl;
    
    t.start();
    for( int i = 0; i < 200000000; i++ )
    { /* empty */
    }
    t.stop();

    t.print( std::cout, 6 );
    
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
    for( int arg = 1; arg < argc; arg++ )
	if( string( argv[arg] ) == "--version" )
	{
	    if( rtt_c4::node() == 0 )
		cout << argv[0] << ": version " << rtt_c4::release() 
		     << endl;
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
	rtt_c4::finalize();
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
