//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Timer.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:56:11 2002
 * \brief  Timer member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iomanip>
#include <c4/config.h>
#include "Timer.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//

//! Constructor
Timer::Timer()
    : begin( 0.0 ),
      end(   0.0 ),
      clock_resolution( wall_clock_resolution() ),
      timer_on(          false ),
      sum_wall( 0.0 ),
      sum_system( 0.0 ),
      sum_user( 0.0 ),
      num_intervals( 0 ),
      component_times_available( set_cta() )
{
    reset();
}

//---------------------------------------------------------------------------//
// Memeber functions
//---------------------------------------------------------------------------//

//! Print out a timing report.
void Timer::print(std::ostream &out, int p) const
{
    using std::setw;
    using std::endl;
    using std::ios;
    
    out.setf(ios::fixed, ios::floatfield);
    out.precision(p);
    out << '\n';
    
    if ( num_intervals > 1 )
	out << "LAST INTERVAL: " << '\n';

    out << setw(20) << "WALL CLOCK TIME: " << wall_clock() << '\n';
    if( component_times_available )
    {
	out << setw(20) << "  USER CPU TIME: " << user_cpu() << '\n';
	out << setw(20) << "SYSTEM CPU TIME: " << system_cpu() << '\n';
    }
    out << '\n';
    
    if ( num_intervals > 1 )
    {
	out << "OVER " << num_intervals << " INTERVALS: " << '\n';
	out << setw(20) << "WALL CLOCK TIME: " << sum_wall_clock() << '\n';
	if( component_times_available )
	{
	    out << setw(20) << "  USER CPU TIME: " << sum_user_cpu() << '\n';
	    out << setw(20) << "SYSTEM CPU TIME: " << sum_system_cpu() << '\n';
	}
	out << '\n';
    }
    out.flush();
}

//---------------------------------------------------------------------------//
// Is this an MPI or Posix timer?
//---------------------------------------------------------------------------//
bool Timer::set_cta()
{
#ifdef C4_SCALAR
    return true;
#else
    return false;
#endif
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                              end of Timer.cc
//---------------------------------------------------------------------------//
