//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Particle.t.hh"
 * \author Thomas M. Evans, Todd J. Urbatsch, Mike Buksas
 * \date   Fri Jan 30 17:04:24 1998
 * \brief  Particle class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Particle.hh"
#include "Global.hh"
#include <iomanip>
#include <algorithm>

namespace rtt_imc 
{

//===========================================================================//
// PARTICLE<MT> FUNCTIONS
//===========================================================================//

//---------------------------------------------------------------------------//
// STATIC MEMBER VARIABLES
//---------------------------------------------------------------------------//

// Parameter minwt_frac is the fractional energy-weight cutoff between
// implicit and analog absorption behavior.
template<class MT> const double Particle<MT>::minwt_frac = 0.01;

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//
// convert a particle event descriptor string into an int

template<class MT>
int Particle<MT>::convert_string_to_descriptor(std::string desc)
{
    // declare return type
    int return_value;

    // born descriptors
    if (desc == "born")
	return_value = BORN;
    else if (desc == "census_born")
	return_value = CENSUS_BORN;
    else if (desc == "boundary_born")
	return_value = BOUNDARY_BORN;
    else if (desc == "vol_emission")
	return_value = VOL_EMISSION;
    else if (desc == "surface_source")
	return_value = SURFACE_SOURCE;
    else if (desc == "unpacked")
	return_value = UNPACKED;

    // collision event descriptors
    else if (desc == "scatter")
	return_value = SCATTER;
    else if (desc == "low_weight")
	return_value = LOW_WEIGHT;
    else if (desc == "eff_scatter")
	return_value = EFF_SCATTER;
    else if (desc == "thom_scatter")
	return_value = THOM_SCATTER;
 
    // streaming descriptors
    else if (desc == "reflection")
	return_value = REFLECTION;
    else if (desc == "stream")
	return_value = STREAM;
    else if (desc == "escape")
	return_value = ESCAPE;
    else if (desc == "cross_boundary")
	return_value = CROSS_BOUNDARY;

    // time and census descriptors
    else if (desc == "census")
	return_value = CENSUS;

    // death
    else if (desc == "killed")
	return_value = KILLED;

    // last else
    else 
	Insist(0,"Invalid string descriptor");

    // return
    return return_value;
}

//---------------------------------------------------------------------------//
// convert an int into a particle event descriptor
   
template<class MT>
std::string Particle<MT>::convert_descriptor_to_string(int index)
{
    std_string desc;

    switch (index) {

	// born descriptors
    case BORN:           
	desc = "born"; 
	break;
    case CENSUS_BORN:    
	desc = "census_born"; 
	break;
    case BOUNDARY_BORN:  
	desc = "boundary_born"; 
	break;
    case VOL_EMISSION:  
	desc = "vol_emission"; 
	break;
    case SURFACE_SOURCE: 
	desc = "surface_source"; 
	break;
    case UNPACKED:
	desc = "unpacked";
	break;

	// collision event descriptors
    case SCATTER:       
	desc = "scatter";
	break;
    case LOW_WEIGHT: 	 
	desc = "low_weight";
	break;
    case EFF_SCATTER:	
	desc = "eff_scatter";
	break;
    case THOM_SCATTER:
	desc = "thom_scatter";
	break;
	
	// streaming descriptors
    case REFLECTION:	
	desc = "reflection";
	break;
    case STREAM:	 
	desc = "stream";
	break;
    case ESCAPE:	
	desc = "escape";
	break;
    case CROSS_BOUNDARY:
	desc = "cross_boundary";
	break;
	
	// time and census descriptors
    case CENSUS:	 
	desc =  "census";
	break;
	
	// death
    case KILLED:	 
	desc =  "killed";
	break;
	
    default:
	Insist(0,"Unrecognized descriptor number");
	break;
    }  

    return desc;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate and return the size of the packed particle.
 * 
 * \param dimension spatial dimension of the problem
 * \param control rtt_rng::Rnd_Control object
 */
template<class MT>
int Particle<MT>::get_packed_particle_size(int dimension,
					   const rtt_rng::Rnd_Control &control)
{
    Require (dimension > 0 && dimension <= 3);

    // determine the size of the packed random number
    int size_rn   = control.get_size();
    Check (size_rn > 0);
    
    // calculate size: int for dimension, int for cell, int for size_rn;
    // dimension doubles for position, 3 doubles for omega, 1 double for ew,
    // 1 double for time_left, 1 double for fraction; size_rn characters for
    // random number state
    int size = 3 * sizeof(int) + (dimension + 6) * sizeof(double) + size_rn;
    
    return size;
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport step for the particle.
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
void Particle<MT>::transport(const MT &mesh, 
			     const Opacity<MT> &xs, 
			     Tally<MT> &tally, 
			     rtt_dsxx::SP<Diagnostic> diagnostic)
{
    Require (alive);

    // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->header();
	diagnostic->print(*this);
    }
  
    // !!! BEGIN TRANSPORT LOOP !!!

    // transport loop, ended when alive = false
    while (alive)
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

	sigma_thomson_scatter = xs.get_sigma_thomson(cell);
	sigma_eff_scatter     = xs.get_sigeffscat(cell);
	sigma_eff_abs         = xs.get_sigeffabs(cell);

	sigma_scatter         = sigma_thomson_scatter + sigma_eff_scatter;

	if (use_analog_absorption())
 	    sigma_analog_abs     = sigma_eff_abs;
	else
	    sigma_analog_abs     = 0.0;

	sigma_collide = sigma_scatter + sigma_analog_abs;

	Check(sigma_collide>=0);

	// accumulate momentum deposition from volume emission particles
	if (descriptor == VOL_EMISSION)
	    tally.accumulate_momentum(cell, -ew, omega);
        
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
	    d_collide = -log(random->ran()) / sigma_collide;

	    prob_thomson_scatter = sigma_thomson_scatter / sigma_collide;
	    prob_scatter         = sigma_scatter         / sigma_collide;
	    prob_abs             = sigma_analog_abs      / sigma_collide;
	}

	Check(d_collide>0);

	// get distance-to-boundary and cell face
	d_boundary = mesh.get_db(r, omega, cell, face);  Check(d_boundary>=0);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * time_left;   Check(d_census);

	// distance until cutoff weight is reached:
	if (sigma_eff_abs == 0 || use_analog_absorption() )
	{ 
	    d_cutoff = rtt_mc::global::huge;
	}
	else
	{
	    d_cutoff = log(fraction/minwt_frac)/sigma_eff_abs;
	}

	Check(d_cutoff>0);


	// detailed diagnostics
	if (diagnostic)
	    if (diagnostic->detail_status())
	    {
		diagnostic->print_dist(d_collide, d_boundary, d_census, cell); 
		diagnostic->print_xs(xs, cell);
	    }


	// determine limiting event
	if      (d_collide < d_boundary  &&  d_collide < d_census   &&
		 d_collide < d_cutoff  )
	{
	    descriptor = COLLISION;
	    d_stream = d_collide;
	}	
	else if (d_boundary < d_collide  &&  d_boundary < d_census  &&
		 d_boundary < d_cutoff )
	{
	    descriptor = BOUNDARY;
	    d_stream = d_boundary;
	}
	else if (d_census < d_collide    &&  d_census < d_boundary  &&  
		 d_census < d_cutoff   )
	{
	    descriptor = CENSUS;
	    d_stream = d_census;
	}
	else if (d_cutoff < d_collide    &&  d_cutoff < d_boundary  &&
		 d_cutoff < d_census   )
	{
	    descriptor = CUTOFF;
	    d_stream = d_cutoff;
	}
	else
	{
	    Insist(0,"D'oh! Transport could not decide limiting event!");
	}


	// Stream the particle, according to its status:
	if (use_analog_absorption())
	{
	    // Light particle (analog) streaming.
	    stream_analog_capture(tally, d_stream);          
	}
	else
	{
	    // Heavy particle (implicit) streaming
	    stream_implicit_capture(xs, tally, d_stream);    
	}


	// Adjust the time remaining till the end of the time step
	time_left -= d_stream / rtt_mc::global::c;


	// Process collisions, boundary crossings, going to census or
	// reaching cutoff events.
	switch (descriptor) 
	{

	case COLLISION:

	    collision_event(mesh, tally, prob_scatter, prob_thomson_scatter, 
			    prob_abs);
	    break;

	case CUTOFF:

	    Check(rtt_mc::global::soft_equiv(fraction, minwt_frac));
	    // Ensure light weight from now on:
	    fraction = minwt_frac * 0.5;  
	    break;

	case BOUNDARY:

	    boundary_event(mesh, tally, face);
	    break;

	case CENSUS:
	    census_event(tally);
	    break;

	}

    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
// PRIVATE TRANSPORT MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
// Process a collision event

template<class MT> 
void Particle<MT>::collision_event(
    const MT& mesh, Tally<MT> &tally, 
    double prob_scatter, double prob_thomson_scatter, double prob_abs){

    double rand_selector = random->ran();
    
    if (rand_selector < prob_scatter) 
    { 
	// Scatter
	if (rand_selector < prob_thomson_scatter)
	{ 
	    // Thomson scatter
	    descriptor = THOM_SCATTER;
	    tally.accum_n_thomscat();
	}
	else
	{ 
	    // Effective scatter
	    descriptor = EFF_SCATTER;
	    tally.accum_n_effscat();
	}
	
	// accumulate momentum from before the scatter
	tally.accumulate_momentum(cell, ew, omega);
	
	// scatter the particle -- update direction cosines
	scatter( mesh );
	
	// accumulate momentum from after the scatter
	tally.accumulate_momentum(cell, -ew, omega);
	
    }
    else if (rand_selector < prob_scatter + prob_abs)
    { // Absorption
	tally.deposit_energy( cell, ew );
	tally.accum_n_killed();
	tally.accum_ew_killed( ew );
	tally.accumulate_momentum(cell, ew, omega);
	descriptor = KILLED; 
	alive=false;
    }
    else
    {
	Insist(0,"D'oh! Transport could not pick a random event!");
    }   
}

//---------------------------------------------------------------------------//
// Process a particle going into census

template<class MT>
void Particle<MT>::census_event(Tally<MT>& tally)
{
    tally.accumulate_cen_info( cell, ew );
    alive = false;
    Check(rtt_mc::global::soft_equiv(time_left, 0.0));
    time_left = 0.0;
}

//---------------------------------------------------------------------------//
// Process a boundary crossing event

template<class MT>
void Particle<MT>::boundary_event(const MT &mesh, Tally<MT>& tally, int face)
{
    tally.accum_n_bndcross();
    alive = surface(mesh, face);
    
    if (descriptor == REFLECTION)
	tally.accum_n_reflections();
    
    if (descriptor == ESCAPE)
    {
	tally.accum_n_escaped();
	tally.accum_ew_escaped(ew);
    }
}

//---------------------------------------------------------------------------//
// perform an isotropic effective scatter 

template<class MT>
void Particle<MT>::scatter(const MT &mesh)
{   
    // calculate theta and phi (isotropic)
    double costheta, phi;
    costheta = 1 - 2 * random->ran();
    phi      = 2 * rtt_mc::global::pi * random->ran();
    
    // get new direction cosines
    mesh.get_Coord().calc_omega(costheta, phi, omega);
}

//---------------------------------------------------------------------------//
// do surface crossings

template<class MT>
bool Particle<MT>::surface(const MT &mesh, int face)
{
    using std::vector;

    // handle particles at a surface

    // status from surface crossing
    bool status;

    // determine the next cell
    int next_cell = mesh.next_cell(cell, face, r);

    // determine descriptor and outcome of this event

    if (next_cell == cell)
    {
	// reflection
	descriptor            = REFLECTION;
	vector<double> normal = mesh.get_normal(cell, face);
	double factor         = rtt_mc::global::dot(omega, normal);
	for (int i = 0; i < mesh.get_Coord().get_sdim(); i++)
	    omega[i] -= 2 * factor * normal[i];
	cell = next_cell;
    }
    else if (next_cell == 0)
    {
	// escape
	descriptor = ESCAPE;
	cell       = next_cell;
    }
    else if (next_cell < 0)
    {
	// domain boundary crossing
	descriptor = CROSS_BOUNDARY;
	cell       = next_cell;
    }
    else 
    {
	// continue streaming
	descriptor = STREAM;
	cell       = next_cell;
    }

    // return outcome of the event
    if (next_cell <= 0)
	status = false;
    else 
	status = true;	    
    return status;
}

//---------------------------------------------------------------------------//
// DIAGNOSTIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
// print out a particle to some stream

template<class MT>
void Particle<MT>::print(std::ostream &output) const
{
    using std::ios;
    using std::setiosflags;
    using std::endl;
    using std::setw;

    // set precisions
    output.precision(3);
    output << setiosflags(ios::fixed);
    
    output << "*** PARTICLE DATA ***" << endl; 
    output << "---------------------" << endl;
    
    // coordinates
    output << setw(20) << setiosflags(ios::right) << "Coordinates: ";
    for (int i = 0; i < r.size(); i++)
	output << setw(12) << r[i] << " ";
    output << endl;
    
    // direction
    output << setw(20) << setiosflags(ios::right) << "Direction: ";
    for (int i = 0; i < omega.size(); i++)
	output << setw(12) << omega[i] << " ";
    output << endl;
    
    // cell
    output << setw(20) << setiosflags(ios::right) << "Cell: " << setw(12) 
	   << cell << endl;
    
    // energy-weight, ew
    output << setw(20) << setiosflags(ios::right) << "Energy-weight: " 
           << setw(12) << ew << endl;
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
// test Particle equality, remember, this simply checks to see if the Sprng
// random number object has the same streamnum, not ID pointer

template<class MT>
bool Particle<MT>::operator==(const Particle<MT> &rhs) const
{
    // check particle data
    if (ew != rhs.ew)
	return false;
    else if (r != rhs.r)
	return false;
    else if (omega != rhs.omega)
	return false;
    else if (cell != rhs.cell)
	return false;
    else if (time_left != rhs.time_left)
	return false;
    else if (fraction != rhs.fraction)
	return false;
    else if (alive != rhs.alive)
	return false;
    else if (descriptor != rhs.descriptor)
	return false;

    if (random->get_num() != rhs.random->get_num())
	return false;

    // if all these things check out then the particles are equal
    return true;
}

//===========================================================================//
// CLASS PARTICLE<MT>::DIAGNOSTIC
//===========================================================================//

template<class MT>
void Particle<MT>::Diagnostic::header() const 
{ 
    output << "*** PARTICLE HISTORY ***" << std::endl; 
    output << "------------------------" << std::endl;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print(const Particle<MT> &particle)  const
{
    using std::ios; 

    // set output precision
    output.precision(3);
    output.setf(ios::scientific, ios::floatfield);

    // print particulars of the particle based on its status
    if (particle.alive == true)
	print_alive(particle);
    else
	print_dead(particle);
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print_alive(const Particle<MT> 
					   &particle) const 
{
    using std::endl;
    using std::ios;
    using std::setiosflags;
    using std::setw;

    // print active particle (alive = true)
    output << " -- Particle is alive -- " << endl;

    // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << convert_descriptor_to_string(particle.descriptor) 
	   << endl;
    
    // coordinates
    output << setw(20) << setiosflags(ios::right) << "Coordinates: ";
    for (int i = 0; i < particle.r.size(); i++)
	output << setw(12) << particle.r[i] << " ";
    output << endl;
    
    // direction
    output << setw(20) << setiosflags(ios::right) << "Direction: ";
    for (int i = 0; i < particle.omega.size(); i++)
	output << setw(12) << particle.omega[i] << " ";
    output << endl;
    
    // cell
    output << setw(20) << setiosflags(ios::right) << "Cell: " << setw(12) 
	   << particle.cell << endl;
    
    // energy-weight, ew
    output << setw(20) << setiosflags(ios::right) << "Energy-weight: " 
           << setw(12) << particle.ew << endl;

    // fraction of original weight
    output << setw(20) << setiosflags(ios::right) << "Fraction: " 
           << setw(12) << particle.fraction << endl;

    // time remaining in this time step
    output << setw(20) << setiosflags(ios::right) << "Time_Left: " 
           << setw(12) << particle.time_left << endl;
    
    output << endl;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print_dead(const Particle<MT> 
					  &particle) const
{
    using std::endl;
    using std::ios;
    using std::setw;
    using std::setiosflags;

    // print dead particle (alive = false)
    output << " -- Particle is dead -- " << endl;

    // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << convert_descriptor_to_string(particle.descriptor) 
	   << endl;
    
    // coordinates
    output << setw(20) << setiosflags(ios::right) << " Last Coordinates: ";
    for (int i = 0; i < particle.r.size(); i++)
	output << setw(12) << particle.r[i] << " ";
    output << endl;
    
    // direction
    output << setw(20) << setiosflags(ios::right) << " Last Direction: ";
    for (int i = 0; i < particle.omega.size(); i++)
	output << setw(12) << particle.omega[i] << " ";
    output << endl;

    // cell
    output << setw(20) << setiosflags(ios::right) << "Last Cell: " 
	   << setw(12) << particle.cell << endl;
    
    // energy-weight, ew
    output << setw(20) << setiosflags(ios::right) << "Last Energy-weight: "
	   << setw(12) << particle.ew << endl;

    // fraction of original weight
    output << setw(20) << setiosflags(ios::right) << "Last Fraction: " 
           << setw(12) << particle.fraction << endl;

    // time remaining in this time step
    output << setw(20) << setiosflags(ios::right) << "Last Time_Left: " 
           << setw(12) << particle.time_left << endl;

    output << endl;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print_dist(double d_scat, double d_bnd, 
					  double d_cen, int cell) const
{
    using std::ios;
    using std::setw;
    using std::endl;

    // do detailed diagnostic print of particle event distances
    output << setw(20) << setiosflags(ios::right) << "Present cell: "
	   << setw(12) << cell << endl;
    output << setw(20) << setiosflags(ios::right) << "Dist-scatter: "
	   << setw(12) << d_scat << endl;
    output << setw(20) << setiosflags(ios::right) << "Dist-boundary: "
	   << setw(12) << d_bnd << endl;   
    output << setw(20) << setiosflags(ios::right) << "Dist-census: "
	   << setw(12) << d_cen << endl;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print_xs(const Opacity<MT> &xs,
					int cell) const
{
    using std::setw;
    using std::endl;
    using std::ios;
    using std::setiosflags;

    // do detailed diagnostic print of particle event cross sections
    output << setw(20) << setiosflags(ios::right) << "Opacity: " 
	   << setw(12) << xs.get_sigma_abs(cell)  << endl;
    output << setw(20) << setiosflags(ios::right) << "Eff. scatter: "
	   << setw(12) << xs.get_sigeffscat(cell) << endl; 
    output << setw(20) << setiosflags(ios::right) << "Eff. absorption: " 
	   << setw(12) << xs.get_sigeffabs(cell)  << endl; 
    output << setw(20) << setiosflags(ios::right) << "Thomson scatter: " 
	   << setw(12) << xs.get_sigma_thomson(cell)  << endl; 
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Particle.t.hh
//---------------------------------------------------------------------------//
