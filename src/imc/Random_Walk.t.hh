//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Random_Walk.t.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 14 18:26:24 2003
 * \brief  Random_Walk member definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Random_Walk_t_hh
#define rtt_imc_Random_Walk_t_hh

#include "Random_Walk.hh"
#include "mc/Constants.hh"
#include <utility>
#include <cmath>

namespace rtt_imc
{

//===========================================================================//
// RANDOM_WALK
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param diff_opacity_in smart pointer to a rtt_imc::Diffusion_Opacity
 * \param min_optical_radius_in number of optical radii required to do a
 * random walk step 
 */
template<class MT>
Random_Walk<MT>::Random_Walk(SP_Mesh              mesh_in,
			     SP_Diffusion_Opacity diff_opacity_in,
			     double               min_optical_radius_in)
    : mesh(mesh_in),
      diff_opacity(diff_opacity_in),
      rw_set(false),
      min_optical_radius(min_optical_radius_in)
{
    Require (diff_opacity);
    Require (mesh);
    Require (mesh->num_cells() == diff_opacity->num_cells());
    Require (min_optical_radius >= 0.0);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief See if we should do Random Walk on this step.
 *
 * Random Walk is valid during a step if:
 * - Ro > min_optical_radius * mean-free-path AND
 * - Ro > distance-to-collision AND
 * - distance-to-census > distance-to-collision
 * .
 *
 * \param cell current cell index
 * \param rw_radius Random Walk sphere radius
 * \param d_collision distance to collision (cm)
 * \param d_census distance to census (shakes)
 */
template<class MT>
bool Random_Walk<MT>::do_a_random_walk(int    cell,
				       double rw_radius_in,
				       double d_collision,
				       double d_census)
{
    Require (cell > 0);
    Require (cell <= diff_opacity->num_cells());
    Require (rw_radius_in >= 0.0);
    Require (rw_set == false);

    // get the rosseland opacity in this cell
    double rosseland = diff_opacity->get_Rosseland_opacity(cell);
    Check (rosseland >= 0.0);

    // calculate the mean free path
    double mfp = 0.0;
    if (rosseland == 0.0)
	mfp = rtt_mc::global::huge;
    else
	mfp = 1.0 / rosseland;

    // set the random walk radius for this step
    rw_radius = rw_radius_in;

    // check to see if random walk is valid

    // rw_radius must be greater than the minimum number of mfps
    if (rw_radius <= min_optical_radius * mfp)
	return false;

    // rw_radius must be greater than the distance to collision
    else if (rw_radius <= d_collision)
	return false;

    // rw_radius is greater than dist-to-collision so rw is possible;
    // however, if d_census is the limiting event then don't do random walk 
    else if (d_census <= d_collision)
	return false;

    // set the rw_set flag to true (random walk is on)
    rw_set = true;

    // return random walk boolean flag
    return true;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Do a random walk.
 * 
 * This function performs a random walk.  The function do_a_random_walk()
 * should already have been called to set the rw_radius (random walk sphere
 * radius).  This is checked in a rtt_dsxx::Require statement.  
 *
 * This function will set a new,
 * - particle direction
 * - particle position
 * - time_left
 * .
 * 
 * \param r particle position; updated by random walk
 * \param omega particle direction; updated by random walk
 * \param time_left time left until end of time step; updated by random walk 
 * \param random particle random number generator
 * \param to_census set to true if random walk particle goes to census
 * \return elapsed time during random walk and distance travelled 
 */
template<class MT>
std::pair<double,double> Random_Walk<MT>::random_walk(sf_double   &r,
						      sf_double   &omega,
						      double      &time_left,
						      int          cell,
						      Rnd_Type     random,
						      bool        &to_census)
{
    Require (rw_set);
    Require (time_left > 0.0);
    Require (rw_radius > 0.0);
    Require (cell > 0);
    Require (cell <= mesh->num_cells());

    // elapsed time during random walk and distance travelled
    pair_double time_radius = std::make_pair(0.0, 0.0);

    // random walk diffusion coefficient
    double D = diff_opacity->get_random_walk_D(cell);
    Check (D > 0.0);

    // check to see if the particle makes it to the surface of the sphere 
    double P_escape = table.get_prob_exit(time_left, D, rw_radius);
    Check (P_escape >= 0.0);
    Check (P_escape <= 1.0);

    // sample to see if particle escapes
    double ran = random.ran();

    // if particle escapes then sample a position on the sphere and update
    // time left
    if (ran <= P_escape)
    {
	// determine the time spent in transport
	time_radius.first = table.get_elapsed_time(D, rw_radius, ran);
	time_left        -= time_radius.first; 
	Check (time_left >= 0.0);

	// the particle sits on the surface of the random walk sphere
	time_radius.second = rw_radius;

	// particle does not goto census
	to_census = false;
    }

    // else determine where the particle is in the sphere
    else
    {
	// we need to determine the radius of the sphere for a particle that
	// goes to census
	
	// first get a new random number
	ran = random.ran();

	// now determine the radius of the sphere the particle lives on
	time_radius.second = table.get_radius(time_left, D, rw_radius, ran);
	Check (time_radius.second >= 0.0);
	Check (time_radius.second <= rw_radius);

	// the particle time_left is zero and the elapsed time is time_left
	time_radius.first = time_left;
	time_left         = 0.0;

	// particle goes to census
	to_census = true;
    }

    // now sample a new particle position on the sphere
    std::pair<sf_double, sf_double> r_and_normal = 
	mesh->sample_pos_on_sphere(cell, r, time_radius.second, random);

    // assign the particle position
    r = r_and_normal.first;
    Check (mesh->in_cell(cell, r));

    // assign the existing particle direction to the normal
    omega = r_and_normal.second;
    
    // sample a new direction from a cosine distribution about the normal
    double costheta = std::sqrt(random.ran());
    double phi      = 2.0 * rtt_mc::global::pi * random.ran();
    mesh->get_Coord().calc_omega(costheta, phi, omega);
    
    // unset random walk
    rw_set = false;

    // return elapsed time
    return time_radius;
}

} // end namespace rtt_imc

#endif                          // rtt_imc_Random_Walk_t_hh

//---------------------------------------------------------------------------//
//                        end of imc/Random_Walk.t.hh
//---------------------------------------------------------------------------//
