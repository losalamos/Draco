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

#include "Timer.hh"
#include <iomanip>

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// PRINTOUT
//---------------------------------------------------------------------------//

void Timer::print(std::ostream &out, int p) const
{
    using std::setw;
    using std::endl;
    using std::cout;
    using std::ios;
    
    out.setf(ios::fixed, ios::floatfield);
    out.precision(p);
    out << endl;
    out << setw(20) << "WALL CLOCK TIME: " << wall_clock() << endl;
    out << setw(20) << "  USER CPU TIME: " << user_cpu() << endl;
    out << setw(20) << "SYSTEM CPU TIME: " << system_cpu() << endl;
    out << endl;
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                              end of Timer.cc
//---------------------------------------------------------------------------//
