//----------------------------------*-C++-*----------------------------------//
// ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the base class time-step advisor.
//---------------------------------------------------------------------------//

#include "timestep/ts_advisor.hh"

#include "c4/global.hh"

#include <iostream>

using std::string; 
using std::cout;
using std::endl;

using namespace rtt_timestep;

ts_advisor::ts_advisor(const string &name_,
		       const usage_flag usage_,
		       const bool active_)

    : name(name_), usage(usage_), active(active_) 
{
// empty
}

ts_advisor::~ts_advisor()
{
// empty
}

void ts_advisor::print(const ts_manager &tsm, const bool controlling) const
{
    if (C4::node() != 0)
	return;
    
    string status = advisor_usable(tsm) ? "true " : "false";
    string space  = "   ";
    string cflag = controlling ? "  ==> " : "      ";
    cout << cflag ;
    cout.width(12);
    cout << get_dt_rec(tsm) << space << 
	status << space << name << endl;
}
    

//---------------------------------------------------------------------------//
//                              end of ts_advisor.cc
//---------------------------------------------------------------------------//
