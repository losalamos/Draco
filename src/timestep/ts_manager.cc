//----------------------------------*-C++-*----------------------------------//
// ts_manager.cc
// John McGhee
// Mon Apr  6 17:22:53 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "timestep/ts_manager.hh"

#include "ds++/Assert.hh"

#include <functional>

#include <stdexcept>

#include <iostream>

using std::list;

typedef ts_advisor TSA;

using std::endl;
using std::cerr;
using std::cout;
using std::ios;

ts_manager::ts_manager()
    :time(0.), dt_new(TSA::small()), dt(TSA::small()), cycle(9999),  
     controlling_advisor("Not Set")
{
    Ensure(invariant_satisfied());
}

ts_manager::~ts_manager()
{
// empty
}

bool ts_manager::invariant_satisfied()
{
    bool ldum =
	0.0 < dt_new &&
	0.0 < dt &&
	controlling_advisor.length() != 0;
    return ldum;
}

double ts_manager::compute_new_timestep(double dt_, int cycle_,
					double time_)
{

    cycle = cycle_;
    dt    = dt_;
    time  = time_;

// Check to be sure that there is at least one usable advisor

    bool found = false;
    for (list< SP<ts_advisor> >::iterator py = advisors.begin(); 
	 py != advisors.end(); py++) 
    {
	if ((**py).advisor_usable(cycle))
	{
	    found = true;
	    break;
	}
    }
    
    if (!found)
    {
	cerr << "  ** Time-Step Manager Warning **" << endl;
	cerr << "  No usable time-step advisors found," << endl;
	cerr << "  defaulting to current time-step" << endl;
	dt_new = dt_;
	controlling_advisor = "Current Time-Step";
	return dt_new;
    }

// Check for one and only one mandatory advisor

    int i = 0;
    for (list< SP<ts_advisor> >::iterator py = advisors.begin(); 
	 py != advisors.end(); py++) 
    {
	if ((**py).advisor_usable(cycle) && (**py).get_usage() == TSA::req)
	{
	    i++;
	    dt_new = (**py).get_dt_rec();
	    controlling_advisor = (**py).get_name();
	}
    } 
    if (i == 1)
    {
	return dt_new;
    }
    else if (i != 0)
    {
	cerr << "  ** Time-Step Manager Warning **" << endl;
	cerr << "  Cycle Number: " << cycle << endl;
	cerr << "  More than one mandatory advisor found," << endl;
	cerr << "  defaulting to last found" << endl;
    }

// Loop over the advisors finding the one that controls

    list< SP<ts_advisor> >::iterator py1=advisors.end();
    list< SP<ts_advisor> >::iterator py2=advisors.end();
    double x1 = TSA::small();
    double x2 = TSA::large();
    for (list< SP<ts_advisor> >::iterator py = advisors.begin(); 
	 py != advisors.end(); py++) 
    {   
	if ((**py).advisor_usable(cycle))
	{
            if ((**py).get_usage() == TSA::min)
	    {
		if ((**py).get_dt_rec() > x1)
		{
		    x1 = (**py).get_dt_rec();
		    py1 = py;
		}
	    }
	    else if ((**py).get_usage() == TSA::max)
	    {
		if ((**py).get_dt_rec() <  x2)
		{
		    x2 = (**py).get_dt_rec();
		    py2 = py;
		}
	    }
	}
    }

    if (py1 == advisors.end() && py2 == advisors.end() )
    {
	cerr << "  ** Time-Step Manager Warning **" << endl;
	cerr << "  Cycle Number: " << cycle << endl;
	cerr << "  No usable time-step advisors found," << endl;
	cerr << "  defaulting to current time-step" << endl;
	dt_new = dt_;
	controlling_advisor = "Current Time-Step";
	return dt_new;
    }
    else if (py1 == advisors.end())
    {
	dt_new = x2;
	controlling_advisor = (**py2).get_name();
    }
    else if (py2 == advisors.end())
    {
	dt_new = x1;
	controlling_advisor = (**py1).get_name();
    }
    else
    {
	if (x1 > x2)
	{
	    dt_new = x1;
	    controlling_advisor = (**py1).get_name();
	    cerr << "  ** Time-Step Manager Warning **" << endl;
	    cerr << "  Cycle Number: " << cycle << endl;
	    cerr << "  No window between min and max advisors," << endl;
	    cerr << "  defaulting to min recommended dt" << endl;
	}
	else
	{
	    dt_new = x2;
	    controlling_advisor = (**py2).get_name();
	}
    }
    
    return dt_new;
}

void ts_manager::add_advisor(const SP<ts_advisor> &new_advisor)
{
    bool not_in_use = true;
    for (list< SP <ts_advisor> >::iterator py = advisors.begin(); 
	 py != advisors.end(); py++) 
    {
	if ((**py).get_name() == (*new_advisor).get_name())
	{
	    not_in_use = false;
	    break;
	}
    }
    if (not_in_use)
    {
	advisors.push_front(new_advisor);
    }
    else
    {
	throw std::runtime_error("Name for requested advisor already in use");
    }
}

void ts_manager::remove_advisor( const SP<ts_advisor> &advisor_to_remove)
{
    for (list< SP<ts_advisor> >::iterator py = advisors.begin(); 
	 py != advisors.end(); py++) 
    {
	if ((**py).get_name() == (*advisor_to_remove).get_name())
	{
	    advisors.erase(py);
	    return;
	}
    }
    throw std::runtime_error("Unable to find requested advisor");
}


struct sptsa_less_than : public std::binary_function< SP<ts_advisor>,
			 SP<ts_advisor>, bool > 
{
    bool operator () (const SP<ts_advisor> &sp_lhs, 
		      const SP<ts_advisor> &sp_rhs) const
    {
	return (*sp_lhs < *sp_rhs);
    }
};

void ts_manager::print_summary() const
{
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(4);
    list < SP<ts_advisor> > temp = advisors;
    temp.sort(sptsa_less_than());
    cout << endl;
    cout << "  *** Time-Step Manager Summary ***" << endl;
    cout << "  Cycle Number         : " << cycle << endl;
    cout << "  Problem Time         : " << time << endl;
    cout << "  Current Time-Step    : " << dt << endl;
    cout << "  Recommended Time-Step: " << dt_new << endl;
    cout << "  Controlling Advisor  : " << controlling_advisor << endl;
    cout << endl;
    cout << "  *** Time-Step Advisor Data Table *** " << endl;
    cout << "  Recommendation   In-Use  Name " << endl;
    for (list< SP<ts_advisor> >::const_iterator py = temp.begin(); 
	 py != temp.end(); py++)
    {
	(**py).print(cycle, (**py).get_name() == controlling_advisor);
    }
    cout << endl;
    cout.setf(0,ios::floatfield);
}

void ts_manager::print_advisors() const
{
    cout << endl;
    cout << "*** Time-Step Manager: Advisor Listing ***" << endl;
    for (list< SP<ts_advisor> >::const_iterator 
	     py = advisors.begin(); py != advisors.end(); py++)
    {
	cout << (**py).get_name() << endl;
    }
    cout << endl; 
}


//---------------------------------------------------------------------------//
//                              end of ts_manager.cc
//---------------------------------------------------------------------------//
