//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_tracker.cc
 * \author Mike Buksas
 * \date   Thu Jun 19 11:33:00 2003
 * \brief  Computes and tallies surface crossings for a collection of surfaces.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Surface_tracker.hh"  
#include "Surface_Tally.hh"

#include <iostream>
#include <cmath>

using std::vector;
using rtt_dsxx::SP;
using rtt_mc::Surface;

namespace rtt_imc
{

Surface_tracker::Surface_tracker(
    const vector<Surface_tracker::SP_Surface>& surfaces_
    )
    : surfaces(surfaces_),
      is_inside(surfaces_.size())
{ /* ... */ }


void Surface_tracker::initialize_status(const vector<double>& position,
					const vector<double>& direction)
{

    int i = 0;
    for (surface_iterator surface = surfaces.begin();
	 surface != surfaces.end();
	 ++surface, ++i)
    {
	is_inside[i] = (*surface)->is_inside(position, direction);
    }

} 

void Surface_tracker::tally_crossings(
    const vector<double>& position,
    const vector<double>& direction,
    double distance, double initial_ew, double sigma,
    Surface_Tally& tally)
{

    Check(position.size()  == 3);
    Check(direction.size() == 3);

    vector<double> final_position(position);
    for (int i=0; i!=3; ++i) final_position[i] += distance * direction[i];

    int i = 0;
    for (surface_iterator surface = surfaces.begin();
	 surface != surfaces.end();
	 ++surface, ++i)
    {

	while (1)
	{
	    bool ends_inside = (*surface)->is_inside(final_position, direction);
	    
	    double crossing_distance =
		(*surface)->distance_to(position, direction, is_inside[i]);

	    Check(distance > 0);

	    if (crossing_distance <= distance || is_inside[i] != ends_inside)
	    {
		
		double crossing_ew = 
		    initial_ew * exp(-sigma * crossing_distance);

		std::cout << "Ray crossed surface " << i 
			  << " with energy weight: " << crossing_ew 
			  << " at distance: " << crossing_distance 
			  << " going " << (is_inside[i] ? "outward" : "inward")
			  << "." << std::endl;

		tally.add_to_tally(direction, is_inside[i], crossing_ew);

		is_inside[i] = !is_inside[i];

	    }

	    if (is_inside[i] == ends_inside) break; // out of while loop

	}

    }
    
}

	     



} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Surface_tracker.cc
//---------------------------------------------------------------------------//
