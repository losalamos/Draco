//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Multigroup_Particle.t.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:31:05 2002
 * \brief  Multigroup_Particle implementation.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Multigroup_Particle_t_hh
#define rtt_imc_Multigroup_Particle_t_hh

#include "Multigroup_Particle.hh"
#include "Random_Walk_Sub_Tally.hh"
#include "ds++/Soft_Equivalence.hh"
#include <utility>

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate and return the size of the packed particle.
 * 
 * \param dimension spatial dimension of the problem
 * \param control rtt_rng::Rnd_Control object
 */
template<class MT>
int Multigroup_Particle<MT>::get_packed_particle_size(
    int                         dimension,
    const rtt_rng::Rnd_Control &control)
{
    Require (dimension > 0 && dimension <= 3);

    // determine the size of the packed random number
    int size_rn = control.get_size();
    Check (size_rn > 0);
    
    // calculate size: int for dimension, int for cell, int for size_rn, int
    // for group_index; dimension doubles for position, 3 doubles for omega,
    // 1 double for ew, 1 double for time_left, 1 double for fraction;
    // size_rn characters for random number state
    int size = 4 * sizeof(int) + (dimension + 6) * sizeof(double) + size_rn;
    
    return size;
}

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor is used to unpack a particle that has been packed with
 * the Particle::pack() function.
 */
template<class MT>
Multigroup_Particle<MT>::Multigroup_Particle(const std::vector<char> &packed)
{
    Require (packed.size() >= 4 * sizeof(int) + 6 * sizeof(double));

    // make an unpacker
    rtt_dsxx::Unpacker u;
    
    // set it
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the spatial dimension of the particle
    int dimension = 0;
    u >> dimension;
    Check (dimension > 0 && dimension <= 3);

    // size the dimension and direction 
    Base::r.resize(dimension);
    Base::omega.resize(3);

    // unpack the position
    for (int i = 0; i < dimension; i++)
	u >> Base::r[i];

    // unpack the rest of the data
    u >> Base::omega[0] >> Base::omega[1] >> Base::omega[2] >> Base::cell 
      >> Base::ew >> Base::time_left >> Base::fraction >> group_index;
    Check (Base::time_left   >= 0.0);
    Check (Base::fraction    >= 0.0);
    Check (Base::cell        >  0);
    Check (Base::ew          >= 0.0);
    Check (group_index       >  0)

    // get the size of the RN state
    int size_rn = 0;
    u >> size_rn;
    Check (size_rn > 0);

    // make a packed rn vector
    std::vector<char> prn(size_rn);

    // unpack the rn state
    for (int i = 0; i < size_rn; i++)
	u >> prn[i];

    // rebuild the rn state
    Base::random = new Rnd_Type(prn);
    Check (Base::random->get_num() >= 0);
    Check (Base::random->get_id());

    // assign the descriptor and status
    Base::descriptor = Base::UNPACKED;
    Base::alive      = true;
    
    Ensure (Base::status());
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the Mulitgroup_Particle.
 *
 * This function transports a particle until it:
 * - leaves the mesh
 * - goes to census
 * - gets killed
 * .
 * The optional arguments specify transport options.  If no optional
 * arguments are given then simple transport is performed as described in
 * \ref mg_trans_method.  If random_walk is active then Random Walk is turned
 * on in thick cells as described in \ref mg_trans_method_walk.
 *
 * \subsection mg_trans_method Straight Multigroup Transport Method
 *
 * Particles undergo implicit absorption until their energy weight drops
 * below 0.01 of their original value. During this time their energy weight
 * decreases exponentialy. Once they attain the cutoff energy weight they are
 * explicitly absorbed and their energy weight remains constant.
 *
 * We unify some of the treatment of the two modes (analog and implicit
 * absorption) by defining sigma_analog_abs to be the actual effective
 * absorption for analog (light) particles and 0.0 for implicit (heavy)
 * particles. The total collision cross section is always
 * sigma_scatter+sigma_eff_abs.
 *
 * To prevent the accumulation of particles with energy weights below the
 * cutoff point, particles stream no further than required to reach cutoff if
 * nothing else happens first. The distance to the cutoff is computed from
 * the current fractional energy weight and the effective absorption in the
 * current cell.
 *
 * \subsection mg_trans_method_walk Multigroup Transport Method with Random
 * Walk
 *
 * Particles undergo implicit absorption until their energy weight drops
 * below 0.01 of their original value. During this time their energy weight
 * decreases exponentialy. Once they attain the cutoff energy weight they are
 * explicitly absorbed and their energy weight remains constant.
 *
 * We unify some of the treatment of the two modes (analog and implicit
 * absorption) by defining sigma_analog_abs to be the actual effective
 * absorption for analog (light) particles and 0.0 for implicit (heavy)
 * particles. The total collision cross section is always
 * sigma_scatter+sigma_eff_abs.
 *
 * To prevent the accumulation of particles with energy weights below the
 * cutoff point, particles stream no further than required to reach cutoff if
 * nothing else happens first. The distance to the cutoff is computed from
 * the current fractional energy weight and the effective absorption in the
 * current cell.
 *
 * Random walk is used as a hybrid diffusion scheme to speed up the
 * transport.  The conditions for Random Walk are:
 * \f[
 * R_{0} > \sigma_{R}
 * \f]
 * \f[
 * R_{0} > d_{\mathrm{collision}} 
 * \f]
 * and 
 * \f[
 * R_{0} < d_{\mathrm{census}}
 * \f]
 * 
 * The Random Walk algorithm is described in Fleck and Canfield, \e J. \e
 * Comp. \e Phys., \b 54, 1984.  
 *
 * \param mesh the mesh object (determined by MT)
 * \param xs the opacity object (multigroup)
 * \param tally the tally object
 * \param random_walk (optional) rtt_dsxx::SP to a Random_Walk object; if the
 * SP is active then random walk is on
 * \param diagnostic (optional) rtt_dsxx::SP to a diagnostic object; if the
 * SP is active then particle tracks will get dumped to standard out.
 */
template<class MT>
void Multigroup_Particle<MT>::transport(const MT          &mesh, 
					const MG_Opacity  &xs, 
					Tally<MT>         &tally,
					SP_Random_Walk     random_walk,
					SP_Surface_tracker surface_tracker,
					SP_Diagnostic      diagnostic)
{
    Require (Base::cell > 0);
    Require (Base::cell <= mesh.num_cells());
    Require (tally.num_cells() == mesh.num_cells());

    // >>> RUN THE TRANSPORT

    // run transport with random walk
    if (random_walk)
	rw_transport(mesh, xs, tally, random_walk, surface_tracker, diagnostic);

    // otherwise run simple transport
    else
	straight_transport(mesh, xs, tally, surface_tracker, diagnostic);

    Ensure (!Base::alive);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a particle into a char stream for communication and
 * persistence. 
 */
template<class MT>
std::vector<char> Multigroup_Particle<MT>::pack() const
{
    Require (Base::omega.size() == 3);

    // make a packer
    rtt_dsxx::Packer p;

    // first pack the random number state
    std::vector<char> prn = Base::random->pack();

    // determine the size of the packed particle: 1 int for cell, + 1 int for
    // size of packed RN state + 1 int for dimension of space + 1 int for
    // group index; dimension + 6 doubles; + size of RN state chars
    int size = 4 * sizeof(int) + (Base::r.size() + 6) * sizeof(double) +
	prn.size();

    // set the packed buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the spatial dimension
    p << static_cast<int>(Base::r.size());
    
    // pack the dimension
    for (int i = 0; i < Base::r.size(); i++)
	p << Base::r[i];
    
    // pack the rest of the data
    p << Base::omega[0] << Base::omega[1] << Base::omega[2] << Base::cell 
      << Base::ew << Base::time_left << Base::fraction << group_index;

    // pack the RN state
    p << static_cast<int>(prn.size());
    for (int i = 0; i < prn.size(); i++)
	p << prn[i];

    Ensure (p.get_ptr() == &packed[0] + size);
    return packed;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the Multigroup_Particle.
 *
 * This transports a particle until it leaves the mesh, goes to census, or
 * gets killed.
 */
template<class MT>
void Multigroup_Particle<MT>::straight_transport(
    const MT           &mesh, 
    const MG_Opacity   &xs, 
    Tally<MT>          &tally, 
    SP_Surface_tracker  surface_tracker,
    SP_Diagnostic       diagnostic)
{
    Require (Base::alive);

    // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->header();
	diagnostic->print(*this);
    }

    // initialize the surface tracker if this feature is requested
    if (surface_tracker) 
	surface_tracker->initialize_status(Base::r, Base::omega);

    // distance to collision, boundary, census and cutoff definitions
    double d_collide, d_boundary, d_census, d_cutoff;
    
    // streaming distance definition
    double d_stream;

    // cell face definition
    int face = 0;

    // event probability definitions
    double prob_thomson_scatter, prob_scatter, prob_abs;

    // intermediate cross section definitions:
    double sigma_thomson_scatter, sigma_eff_scatter, sigma_eff_abs, 
	sigma_scatter, sigma_analog_abs, sigma_collide;
  
    // !!! BEGIN TRANSPORT LOOP !!!

    // transport loop, ended when alive = false
    while (Base::alive)
    {
	// initialize cell face
	face = 0; 

	// get cross sections
	sigma_thomson_scatter = xs.get_sigma_thomson(Base::cell, group_index);
	sigma_eff_scatter     = xs.get_sigeffscat(Base::cell, group_index);
	sigma_eff_abs         = xs.get_sigeffabs(Base::cell, group_index);

	// calculate total scattering cross section
	sigma_scatter = sigma_thomson_scatter + sigma_eff_scatter;

	// determine analog absorption cross section
	if (Base::use_analog_absorption())
 	    sigma_analog_abs = sigma_eff_abs;
	else
	    sigma_analog_abs = 0.0;

	// determine total collision cross section
	sigma_collide = sigma_scatter + sigma_analog_abs;
	Check(sigma_collide >= 0);

	// accumulate momentum deposition from volume emission particles
	if (Base::descriptor == Base::VOL_EMISSION)
	    tally.accumulate_momentum(Base::cell, -Base::ew, Base::omega);
        
	// sample distance-to-scatter/absorption (effective scatter or hardball)
	if (sigma_collide == 0 ) 
	{
	    d_collide = rtt_mc::global::huge;    

	    prob_thomson_scatter = 0.0;
	    prob_scatter         = 0.0;
	    prob_abs             = 0.0;
	}
	else 
	{
	    d_collide = -std::log(Base::random->ran()) / sigma_collide;

	    prob_thomson_scatter = sigma_thomson_scatter / sigma_collide;
	    prob_scatter         = sigma_scatter         / sigma_collide;
	    prob_abs             = sigma_analog_abs      / sigma_collide;
	}

	Check (d_collide > 0);

	// get distance-to-boundary and cell face
	d_boundary = mesh.get_db(Base::r, Base::omega, Base::cell, face); 
	Check (d_boundary >= 0);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * Base::time_left;  
	Check (d_census);

	// distance until cutoff weight is reached:
	if (sigma_eff_abs == 0 || Base::use_analog_absorption() )
	{ 
	    d_cutoff = rtt_mc::global::huge;
	}
	else
	{
	    d_cutoff = std::log(Base::fraction / Base::minwt_frac) / 
		sigma_eff_abs;
	}

	Check(d_cutoff > 0);

	// detailed diagnostics
	if (diagnostic)
	    if (diagnostic->detail_status())
	    {
		diagnostic->print_dist(d_collide, d_boundary, d_census,
				       Base::cell); 
		diagnostic->print_xs(xs, Base::cell, group_index);
	    }

	// determine limiting event
	if      (d_collide < d_boundary  &&  d_collide < d_census   &&
		 d_collide < d_cutoff  )
	{
	    Base::descriptor = Base::COLLISION;
	    d_stream         = d_collide;
	}	
	else if (d_boundary < d_collide  &&  d_boundary < d_census  &&
		 d_boundary < d_cutoff )
	{
	    Base::descriptor = Base::BOUNDARY;
	    d_stream         = d_boundary;
	}
	else if (d_census < d_collide    &&  d_census < d_boundary  &&  
		 d_census < d_cutoff   )
	{
	    Base::descriptor = Base::CENSUS;
	    d_stream         = d_census;
	}
	else if (d_cutoff < d_collide    &&  d_cutoff < d_boundary  &&
		 d_cutoff < d_census   )
	{
	    Base::descriptor = Base::CUTOFF;
	    d_stream         = d_cutoff;
	}
	else
	{
	    throw rtt_dsxx::assertion(
		"Transport could not decide limiting event!");
	}

	// Stream the particle, according to its status:
	Base::stream_and_capture(tally, surface_tracker, sigma_eff_abs,
				 d_stream, group_index);

	// Process collisions, boundary crossings, going to census or
	// reaching cutoff events.
	switch (Base::descriptor) 
	{

	case Base::COLLISION:

	    // process the collision event
	    collision_event(mesh, tally, xs, prob_scatter, prob_thomson_scatter,
			    prob_abs);
	    break;

	case Base::CUTOFF:

	    Check(rtt_mc::global::soft_equiv(Base::fraction, 
					     Base::minwt_frac));

	    // Ensure light weight from now on:
	    Base::fraction = Base::minwt_frac * 0.5;  
	    break;

	case Base::BOUNDARY:

	    // process a boundary event
	    Base::boundary_event(mesh, tally, face);
	    break;

	case Base::CENSUS:
	    
	    // process a census event
	    Base::census_event(tally);
	    break;

	default:

	    // throw an assertion
	    throw 
		rtt_dsxx::assertion("Undefined event in Multigroup_Particle.");
	}

	// initialize diagnostics
	if (diagnostic)
	    diagnostic->print(*this);
    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the Multigroup_Particle with Random Walk.
 *
 * This transports a particle until it leaves the mesh, goes to census, or
 * gets killed. Random Walk is used in thick, diffusive cells.
 */
template<class MT>
void Multigroup_Particle<MT>::rw_transport(
    const MT          &mesh, 
    const MG_Opacity  &xs, 
    Tally<MT>         &tally, 
    SP_Random_Walk     random_walk,
    SP_Surface_tracker surface_tracker,
    SP_Diagnostic      diagnostic)
{
    Require (Base::alive);
    Require (tally.get_RW_Sub_Tally());

    // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->header();
	diagnostic->print(*this);
    }

    // initialize the surface tracker if this feature is requested
    if (surface_tracker) 
	surface_tracker->initialize_status(Base::r, Base::omega);

    // distance to collision, boundary, census and cutoff definitions
    double d_collide, d_boundary, d_census, d_cutoff;
    
    // streaming distance definition
    double d_stream;

    // random walk sphere radius
    double rw_radius;

    // cell face definition
    int face = 0;

    // event probability definitions
    double prob_thomson_scatter, prob_scatter, prob_abs;

    // intermediate cross section definitions:
    double sigma_thomson_scatter, sigma_eff_scatter, sigma_eff_abs, 
	sigma_scatter, sigma_analog_abs, sigma_collide;

    // are we doing random walk on this step
    bool do_a_random_walk;
  
    // !!! BEGIN TRANSPORT LOOP !!!

    // transport loop, ended when alive = false
    while (Base::alive)
    {
	// initialize cell face
	face = 0; 

	// get cross sections
	sigma_thomson_scatter = xs.get_sigma_thomson(Base::cell, group_index);
	sigma_eff_scatter     = xs.get_sigeffscat(Base::cell, group_index);
	sigma_eff_abs         = xs.get_sigeffabs(Base::cell, group_index);

	// calculate total scattering cross section
	sigma_scatter = sigma_thomson_scatter + sigma_eff_scatter;

	// determine analog absorption cross section
	if (Base::use_analog_absorption())
 	    sigma_analog_abs = sigma_eff_abs;
	else
	    sigma_analog_abs = 0.0;

	// determine total collision cross section
	sigma_collide = sigma_scatter + sigma_analog_abs;
	Check (sigma_collide >= 0);

	// accumulate momentum deposition from volume emission particles
	if (Base::descriptor == Base::VOL_EMISSION)
	    tally.accumulate_momentum(Base::cell, -Base::ew, Base::omega);
        
	// sample distance-to-scatter/absorption (effective scatter or hardball)
	if (sigma_collide == 0 ) 
	{
	    d_collide = rtt_mc::global::huge;    

	    prob_thomson_scatter = 0.0;
	    prob_scatter         = 0.0;
	    prob_abs             = 0.0;
	}
	else 
	{
	    d_collide = -std::log(Base::random->ran()) / sigma_collide;

	    prob_thomson_scatter = sigma_thomson_scatter / sigma_collide;
	    prob_scatter         = sigma_scatter         / sigma_collide;
	    prob_abs             = sigma_analog_abs      / sigma_collide;
	}
	Check (d_collide > 0);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * Base::time_left;  
	Check (d_census);

	// distance until cutoff weight is reached:
	if (sigma_eff_abs == 0 || Base::use_analog_absorption() )
	{ 
	    d_cutoff = rtt_mc::global::huge;
	}
	else
	{
	    d_cutoff = std::log(Base::fraction / Base::minwt_frac) / 
		sigma_eff_abs;
	}
	Check (d_cutoff > 0);

	// check to see if we should do a random walk; check on initial
	// conditions first
	if (Base::descriptor == Base::VOL_EMISSION || 
	    Base::descriptor == Base::EFF_SCATTER)
	{
	    // get the size of the random walk sphere
	    rw_radius = mesh.get_random_walk_sphere_radius(Base::r, Base::cell);
	    Check (rw_radius >= 0.0);

	    // tally rw sphere radius
	    tally.get_RW_Sub_Tally()->accum_sphere_radii(rw_radius);

	    // check to see if the random walk conditions are valid
	    do_a_random_walk = random_walk->do_a_random_walk(
		Base::cell, rw_radius, d_collide, d_census, xs);

	    if (surface_tracker && surface_tracker->surface_in_cell(Base::cell))
		do_a_random_walk = false;
	}
	else
	{
	    // initial conditions not met, don't do a random walk
	    do_a_random_walk = false;
	}
	 
	// calculate distance to boundary if we are not doing random walk
	if (!do_a_random_walk)
	{
	    // get distance-to-boundary and cell face
	    d_boundary = mesh.get_db(Base::r, Base::omega, Base::cell, face);  
	    Check (d_boundary >= 0);
	}
	else   
	{
	    // we don't need distance to boundary so set it to zero
	    d_boundary = 0.0;
	}

	// detailed diagnostics
	if (diagnostic)
	    if (diagnostic->detail_status())
	    {
		diagnostic->print_dist(d_collide, d_boundary, d_census,
				       Base::cell); 
		diagnostic->print_xs(xs, Base::cell, group_index);
	    }

	// determine limiting event
	if      (do_a_random_walk)
	{
	    Base::descriptor = Base::RANDOM_WALK;
	    d_stream         = 0.0;
	}
	else if (d_collide < d_boundary  &&  d_collide < d_census   &&
		 d_collide < d_cutoff  )
	{
	    Base::descriptor = Base::COLLISION;
	    d_stream         = d_collide;
	}	
	else if (d_boundary < d_collide  &&  d_boundary < d_census  &&
		 d_boundary < d_cutoff )
	{
	    Base::descriptor = Base::BOUNDARY;
	    d_stream         = d_boundary;
	}
	else if (d_census < d_collide    &&  d_census < d_boundary  &&  
		 d_census < d_cutoff   )
	{
	    Base::descriptor = Base::CENSUS;
	    d_stream         = d_census;
	}
	else if (d_cutoff < d_collide    &&  d_cutoff < d_boundary  &&
		 d_cutoff < d_census   )
	{
	    Base::descriptor = Base::CUTOFF;
	    d_stream         = d_cutoff;
	}
	else
	{
	    throw rtt_dsxx::assertion(
		"Transport could not decide limiting event!");
	}

	// stream the particle: transport or random walk
	if (do_a_random_walk)
	{
	    // boolean to see if particle goes to census
	    bool to_census;

	    // do a random walk
	    std::pair<double,double> time_radius = random_walk->random_walk(
		Base::r, Base::omega, Base::time_left, Base::cell, 
		*Base::random, to_census);

	    // sample a new energy group (add 1 since return value is in
	    // [0,G-1])
	    group_index = rtt_mc::sampler::sample_bin_from_discrete_cdf(
		*Base::random, xs.get_emission_group_cdf(Base::cell)) + 1;
	    Check (group_index > 0);
	    Check (group_index <= xs.get_Frequency()->get_num_groups());
	    
	    // tally optical length of random walk
	    tally.get_RW_Sub_Tally()->accum_step_length(
		time_radius.second * sigma_collide);

	    // set descriptor if particle goes to census
	    if (to_census)
		Base::descriptor = Base::CENSUS;

	    // calculate the plank cross section for adjusting weight of
	    // particle after random walk
	    double planck = xs.get_integrated_sigma_times_Planck(Base::cell) / 
		xs.get_integrated_norm_Planck(Base::cell);
	    Check (planck >= 0.0);

	    // process the random walk absorption
	    Base::random_walk_event(time_radius.first, time_radius.second,
				    tally, planck, xs.get_fleck(Base::cell));
	}
	else
	{
	    // Stream the particle, according to its status:
	    Base::stream_and_capture(tally, surface_tracker, sigma_eff_abs,
				     d_stream, group_index);
	}

	// Process collisions, boundary crossings, going to census or
	// reaching cutoff events.
	switch (Base::descriptor) 
	{

	case Base::KILLED:
	    
	    Check (!Base::alive);
	    break;

	case Base::RANDOM_WALK:

	    // we have already processed the event
	    break;

	case Base::COLLISION:

	    // process the collision event
	    collision_event(mesh, tally, xs, prob_scatter, prob_thomson_scatter,
			    prob_abs);
	    break;

	case Base::CUTOFF:

	    Check(rtt_mc::global::soft_equiv(Base::fraction,
					     Base::minwt_frac));

	    // Ensure light weight from now on:
	    Base::fraction = Base::minwt_frac * 0.5;  
	    break;

	case Base::BOUNDARY:

	    // process a boundary event
	    Base::boundary_event(mesh, tally, face);
	    break;

	case Base::CENSUS:
	    
	    // process a census event
	    Base::census_event(tally);
	    break;

	default:

	    // throw an assertion
	    throw 
		rtt_dsxx::assertion("Undefined event in Multigroup_Particle.");
	}

	// initialize diagnostics
	if (diagnostic)
	    diagnostic->print(*this);
    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a collision.
 */
template<class MT> 
void Multigroup_Particle<MT>::collision_event(
    const MT         &mesh, 
    Tally<MT>        &tally, 
    const MG_Opacity &opacity,
    double            prob_scatter, 
    double            prob_thomson_scatter, 
    double            prob_abs)
{
    Check (mesh.num_cells() == tally.num_cells());

    // get a random number
    double rand_selector = Base::random->ran();
    
    if (rand_selector < prob_scatter) 
    { 
	// accumulate momentum from before the scatter
	tally.accumulate_momentum(Base::cell, Base::ew, Base::omega);
	
	// Scatter
	if (rand_selector < prob_thomson_scatter)
	{ 
	    // Thomson scatter
	    Base::descriptor = Base::THOM_SCATTER;
	    tally.accum_n_thomscat();
	
	    // scatter the particle -- update direction cosines
	    Base::scatter(mesh);
	}
	else
	{ 
	    // Effective scatter
	    Base::descriptor = Base::EFF_SCATTER;
	    tally.accum_n_effscat();

	    // scatter the particle -- update direction cosines
	    effective_scatter(mesh, opacity);
	}
	
	// accumulate momentum from after the scatter
	tally.accumulate_momentum(Base::cell, -Base::ew, Base::omega);
    }
    else if (rand_selector < prob_scatter + prob_abs)
    { 
	// Absorption
	
	// tally absorption data
	tally.deposit_energy(Base::cell, Base::ew);
	tally.accum_n_killed();
	tally.accum_ew_killed(Base::ew);
	tally.accumulate_momentum(Base::cell, Base::ew, Base::omega);

	// set the descriptor and particle status
	Base::descriptor = Base::KILLED; 
	Base::alive      = false;
    }
    else
    {
	Insist(0,"D'oh! Transport could not pick a random event!");
    }   
}

//---------------------------------------------------------------------------//
// DIAGNOSTIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out a Multigroup_Particle to some stream.
 */
template<class MT>
void Multigroup_Particle<MT>::print(std::ostream &out) const
{
    using std::ios;
    using std::setiosflags;
    using std::endl;
    using std::setw;

    // call the base class print
    Base::print(out);

    // add the group index to the list
    out << setw(20) << setiosflags(ios::right) << "Group Index: " 
	<< setw(12) << group_index << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test equality of two Multigroup_Particle objects.
 *
 * Test Multigroup_Particle equality.  With regard to the random number state
 * of the particle, this simply checks to see if the Sprng random number
 * object has the same streamnum, not ID pointer
 */
template<class MT>
bool Multigroup_Particle<MT>::operator==(const Multigroup_Particle<MT> &rhs)
    const
{
    // check particle data
    if (Base::ew != rhs.Base::ew)
	return false;
    else if (Base::r != rhs.Base::r)
	return false;
    else if (Base::omega != rhs.Base::omega)
	return false;
    else if (Base::cell != rhs.Base::cell)
	return false;
    else if (Base::time_left != rhs.Base::time_left)
	return false;
    else if (Base::fraction != rhs.Base::fraction)
	return false;
    else if (Base::alive != rhs.Base::alive)
	return false;
    else if (Base::descriptor != rhs.Base::descriptor)
	return false;
    else if (group_index != rhs.group_index)
	return false;

    if (Base::random->get_num() != rhs.Base::random->get_num())
	return false;

    // if all these things check out then the particles are equal
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test inequality of two Gray_Particle objects.
 */
template<class MT>
bool Multigroup_Particle<MT>::operator!=(const Multigroup_Particle<MT> &rhs) 
    const
{
    return !(*this == rhs);
}

//===========================================================================//
// MULTIGROUP_PARTICLE::DIAGNOSTIC FUNCTIONS
//===========================================================================//

template<class MT>
void Multigroup_Particle<MT>::Diagnostic::print_xs(
    const MG_Opacity &xs,
    int               cell_in,
    int               group_index) const
{
    using std::setw;
    using std::endl;
    using std::ios;
    using std::setiosflags;

    // do detailed diagnostic print of particle event cross sections
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right)
			     << "Opacity: " << setw(12) 
			     << xs.get_sigma_abs(cell_in, group_index)
			     << endl;
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right)
			     << "Eff. scatter: " << setw(12) 
			     << xs.get_sigeffscat(cell_in, group_index)
			     << endl; 
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Eff. absorption: " << setw(12)
			     << xs.get_sigeffabs(cell_in, group_index)
			     << endl; 
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right)
			     << "Thomson scatter: " << setw(12)
			     << xs.get_sigma_thomson(cell_in, group_index)  
			     << endl; 
}

} // end namespace rtt_imc

#endif                          // rtt_imc_Multigroup_Particle_t_hh

//---------------------------------------------------------------------------//
//                        end of imc/Multigroup_Particle.t.hh
//---------------------------------------------------------------------------//
