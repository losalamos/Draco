//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_tracker.cc
 * \author Mike Buksas
 * \date   Thu Jun 19 11:33:00 2003
 * \brief  Implementation file for Surface_tracker
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Surface_tracker.hh"  
#include "Surface_Sub_Tally.hh"

#include <iostream>
#include <cmath>
#include <functional>

using namespace std;
using rtt_dsxx::SP;
using rtt_mc::Surface;

namespace rtt_imc
{

Surface_tracker::Surface_tracker(
    const vector<Surface_tracker::SP_Surface>& surfaces_)
    : surface_list(surfaces_),
      is_inside(surfaces_.size()),
      tally_indices(surfaces_.size())
{ 

    int index = 1;
    for (vector<int>::iterator tally_index = tally_indices.begin();
	 tally_index != tally_indices.end(); 
	 ++tally_index)
    {
	*tally_index = index++;
    }
	

}

Surface_tracker::Surface_tracker(
    const vector<Surface_tracker::SP_Surface>& surfaces_,
    const vector<int>& tally_indices_)
    : surface_list(surfaces_),
      is_inside(surfaces_.size()),
      tally_indices(tally_indices_)
{

    // Make sure index list is the same size:
    Require( tally_indices.size() == surface_list.size() );

    // Make sure tally indices are positive:
    Require (
	find_if( tally_indices.begin(), tally_indices.end(),
		 bind2nd(less_equal<int>(), 0) ) == 
	tally_indices.end() 
	) ;
    
    // Make sure tally indices are increasing:
    Require (
	adjacent_find( tally_indices.begin(), tally_indices.end(),
		       greater_equal<int>()) == 
	tally_indices.end() 
	);
    
}


//---------------------------------------------------------------------------//
/*! 
 * \brief Resets the offical inside/outside status variables. Use when
 * tracking for a new particle.
 * 
 * \param position Position of the Particle.
 * \param direction Direction of the Particle.
 * \return void
 *
 * This function calls is_inside for each surface with the provided position
 * and direction. It stores the result of each call in the internal boolean
 * array. 
 */

void Surface_tracker::initialize_status(const std::vector<double>& position,
					const std::vector<double>& direction)
{

    vector<bool>::iterator status = is_inside.begin();
    for (surface_iterator surface = surface_list.begin();
	 surface != surface_list.end();
	 ++surface, ++status)
    {
	*status = (*surface)->is_inside(position, direction);
    }

} 


//---------------------------------------------------------------------------//
/*! 
 * \brief Checks for and tallies surface crossings during an implicit
 * streaming step.
 *
 * This function takes an initial position and direction and a streaming
 * distance and determines whrether any of the collection of surfaces was
 * crossed. For each crossing that takes place, the energy weight at the
 * point of crossing is computed from the initial weight \f$ w\f$ and the
 * coefficient of attenuation \f$\sigma\f$ via \f$ w_c = e^{-\sigma d_c}\f$
 * where \f$ d_c\f$ is the distance to the crossing.
 * 
 * \param position Position of the Particle
 * \param direction Direction of the Particle
 * \param distance distance of travel
 * \param initial_ew energy weight at beginning of step
 * \param sigma coefficient of exponential attenuation with distance traveled.
 * \return void
 */

void Surface_tracker::tally_crossings_implicit_abs(
    const std::vector<double>& position,
    const std::vector<double>& direction,
    double distance, double initial_ew, double sigma,
    Surface_Sub_Tally& tally)
{

    Check(position.size()  == 3);
    Check(direction.size() == 3);

    vector<double> final_position(position);
    for (int i=0; i!=3; ++i) final_position[i] += distance * direction[i];

    int surface_index = 0;
    for (surface_iterator surface = surface_list.begin();
	 surface != surface_list.end();
	 ++surface, ++surface_index)
    {

	bool ends_inside = (*surface)->is_inside(final_position, direction);
	    
	do // while is_inside != ends_inside
	{
	    double crossing_distance =
		(*surface)->distance_to(position, direction, 
					is_inside[surface_index]);

	    Check(distance > 0);

	    if (crossing_distance <= distance || is_inside[surface_index] != 
		ends_inside)
	    {
		
		double crossing_ew = 
		    initial_ew * exp(-sigma * crossing_distance);
		
		tally.add_to_tally(tally_indices[surface_index], 
				   direction, is_inside[surface_index], 
				   crossing_ew);

		is_inside[surface_index] = !is_inside[surface_index];

	    }

	} while (is_inside[surface_index] != ends_inside);

    }
    
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Checks for and tallies surface crossings during an analog streaming
 * step.
 * 
 * This function takes an initial position and direction and a streaming
 * distance and determines whrether any of the collection of surfaces was
 * crossed. Each crossing is tallied at the given energy weight.
 *
 * \param position Position of the Particle
 * \param direction Direction of the Particle
 * \param distance distance of travel
 * \param ew energy weight
 * \return void
 */

void Surface_tracker::tally_crossings_analog_abs(
    const std::vector<double>& position,
    const std::vector<double>& direction,
    double distance, double ew, 
    Surface_Sub_Tally& tally)
{

    Check(position.size()  == 3);
    Check(direction.size() == 3);

    vector<double> final_position(position);
    for (int i=0; i!=3; ++i) final_position[i] += distance * direction[i];

    int surface_index = 0;
    for (surface_iterator surface = surface_list.begin();
	 surface != surface_list.end();
	 ++surface, ++surface_index)
    {

	bool ends_inside = (*surface)->is_inside(final_position, direction);
	    
	do // is_inside != ends_inside
	{
	    double crossing_distance =
		(*surface)->distance_to(position, direction, 
					is_inside[surface_index]);

	    Check(distance > 0);

	    if (crossing_distance <= distance || 
		is_inside[surface_index] != ends_inside)
	    {
		
		tally.add_to_tally(tally_indices[surface_index], 
				   direction, is_inside[surface_index], ew);

		is_inside[surface_index] = !is_inside[surface_index];

	    }

	} while (is_inside[surface_index] != ends_inside);

    }
    
}

	     



} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Surface_tracker.cc
//---------------------------------------------------------------------------//
