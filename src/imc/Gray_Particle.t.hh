//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Gray_Particle.t.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:30:27 2002
 * \brief  Gray_Particle class implementation.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Gray_Particle.hh"

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
int Gray_Particle<MT>::get_packed_particle_size(
    int                         dimension,
    const rtt_rng::Rnd_Control &control)
{
    Require (dimension > 0 && dimension <= 3);

    // determine the size of the packed random number
    int size_rn = control.get_size();
    Check (size_rn > 0);
    
    // calculate size: int for dimension, int for cell, int for size_rn;
    // dimension doubles for position, 3 doubles for omega, 1 double for ew,
    // 1 double for time_left, 1 double for fraction; size_rn characters for
    // random number state
    int size = 3 * sizeof(int) + (dimension + 6) * sizeof(double) + size_rn;
    
    return size;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the Gray_Particle.
 *
 * This transports a particle until it leaves the mesh, goes to census, or
 * gets killed.
 *
 * \p Transport Method
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
 */
template<class MT>
void Gray_Particle<MT>::transport(
    const MT                         &mesh, 
    const Opacity<MT,Gray_Frequency> &xs, 
    Tally<MT>                        &tally, 
    SP_Diagnostic                     diagnostic)
{
    Require (Base::alive);

    // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->header();
	diagnostic->print(*this);
    }
  
    // !!! BEGIN TRANSPORT LOOP !!!

    // transport loop, ended when alive = false
    while (Base::alive)
    {
	// distance to collision, boundary, census and cutoff defnintions
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

	sigma_thomson_scatter = xs.get_sigma_thomson(Base::cell);
	sigma_eff_scatter     = xs.get_sigeffscat(Base::cell);
	sigma_eff_abs         = xs.get_sigeffabs(Base::cell);

	sigma_scatter         = sigma_thomson_scatter + sigma_eff_scatter;

	if (Base::use_analog_absorption())
 	    sigma_analog_abs     = sigma_eff_abs;
	else
	    sigma_analog_abs     = 0.0;

	sigma_collide = sigma_scatter + sigma_analog_abs;

	Check(sigma_collide>=0);

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

	Check(d_collide > 0);

	// get distance-to-boundary and cell face
	d_boundary = mesh.get_db(Base::r, Base::omega, Base::cell, face);  
	Check(d_boundary>=0);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * Base::time_left;   
	Check(d_census);

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
		diagnostic->print_xs(xs, Base:: cell);
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
	    Insist(0,"D'oh! Transport could not decide limiting event!");
	}


	// Stream the particle, according to its status:
	if (Base::use_analog_absorption())
	{
	    // Light particle (analog) streaming.
	    Base::stream_analog_capture(tally, d_stream);          
	}
	else
	{
	    // Heavy particle (implicit) streaming
	    stream_implicit_capture(xs, tally, d_stream);    
	}


	// Adjust the time remaining till the end of the time step
	Base::time_left -= d_stream / rtt_mc::global::c;


	// Process collisions, boundary crossings, going to census or
	// reaching cutoff events.
	switch (Base::descriptor) 
	{

	case Base::COLLISION:

	    // process the collision
	    collision_event(mesh, tally, prob_scatter, prob_thomson_scatter, 
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
	    throw rtt_dsxx::assertion("Undefined event in Gray_Particle.");
	}
    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
// DIAGNOSTIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Test equality of two Gray_Particle objects.
 *
 * Test Gray_Particle equality. With regard to the random number state of the
 * particle, this simply checks to see if the Sprng random number object has
 * the same streamnum, not ID pointer.
 */
template<class MT>
bool Gray_Particle<MT>::operator==(const Gray_Particle<MT> &rhs) const
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
bool Gray_Particle<MT>::operator!=(const Gray_Particle<MT> &rhs) const
{
    return !(*this == rhs);
}

//===========================================================================//
// GRAY_PARTICLE::DIAGNOSTIC FUNCTIONS
//===========================================================================//

template<class MT>
void Gray_Particle<MT>::Diagnostic::print_xs(
    const Opacity<MT,Gray_Frequency> &xs,
    int                               cell_in) const
{
    using std::setw;
    using std::endl;
    using std::ios;
    using std::setiosflags;

    // do detailed diagnostic print of particle event cross sections
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Opacity: " << setw(12) 
			     << xs.get_sigma_abs(cell_in) << endl;
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Eff. scatter: " << setw(12) 
			     << xs.get_sigeffscat(cell_in) << endl; 
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Eff. absorption: " << setw(12) 
			     << xs.get_sigeffabs(cell_in) << endl; 
    Base::Diagnostic::output << setw(20) << setiosflags(ios::right) 
			     << "Thomson scatter: " << setw(12) 
			     << xs.get_sigma_thomson(cell_in) << endl; 
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Gray_Particle.t.hh
//---------------------------------------------------------------------------//
