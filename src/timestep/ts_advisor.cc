//----------------------------------*-C++-*----------------------------------//
// ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "timestep/ts_advisor.hh"

#include <iostream>

using std::string; 
using std::cout;
using std::endl;

ts_advisor::ts_advisor(const string &name_,
		       const usage_flag usage_,
		       const bool active_)

    : name(name_), usage(usage_), active(active_), 
      cycle_at_last_update(-989898), dt_rec(large())
{
// empty
}

ts_advisor::~ts_advisor()
{
// empty
}

void ts_advisor::print(const int cycle_, const bool controlling) const
{

    string status = advisor_usable(cycle_) ? "true " : "false";
    string space  = "   ";
    string cflag = controlling ? "  ==> " : "      ";
    cout << cflag << dt_rec << space << 
	status << space << name << endl;
}
    

//---------------------------------------------------------------------------//
//                              end of ts_advisor.cc
//---------------------------------------------------------------------------//
