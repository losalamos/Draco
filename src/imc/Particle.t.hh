//----------------------------------*-C++-*----------------------------------//
// Particle.t.hh
// Thomas M. Evans
// Fri Jan 30 17:04:24 1998
//---------------------------------------------------------------------------//
// @> Particle class implementation file
//---------------------------------------------------------------------------//

#include "Particle.hh"
#include "Global.hh"
#include <iomanip>

namespace rtt_imc 
{

// STL functions used
using std::endl;
using std::setw;
using std::ios;
using std::setiosflags;
using std::ostream;
using std::vector;

// draco stuff
using rtt_dsxx::SP;

// services from IMC::Global namespace
using rtt_mc::global::pi;
using rtt_mc::global::c;
using rtt_mc::global::dot;

//===========================================================================//
// class Particle<MT>
//===========================================================================//

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

// Particle uses default copy constructors and assignment operators

//---------------------------------------------------------------------------//
// public transport member functions
//---------------------------------------------------------------------------//
// transport a particle

template<class MT>
void Particle<MT>::transport(const MT &mesh, const Opacity<MT> &xs, 
			     Tally<MT> &tally, SP<Diagnostic> diagnostic)
{
    // transport particle through mesh using regular IMC transport
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
	// dist-to-scatter, dist-to-boundary, and dist-to-census definitions
        double d_scatter, d_boundary, d_census;
	double dist_stream;
	// cell face definition
        int face = 0;

	// accumulate momentum deposition from volume emission particles
	if (descriptor == "vol_emission")
	    tally.accumulate_momentum(cell, -ew, omega);
        
	// sample distance-to-scatter (effective scatter or hardball)
	d_scatter = -log(random.ran()) / 
	    (xs.get_sigeffscat(cell) + xs.get_sigma_thomson(cell));

	// get distance-to-boundary and cell face
        d_boundary = mesh.get_db(r, omega, cell, face);

	// distance to census (end of time step)
	d_census = rtt_mc::global::c * time_left;

	// detailed diagnostics
	if (diagnostic)
	    if (diagnostic->detail_status())
	    {
		diagnostic->print_dist(d_scatter, d_boundary, d_census,
				       cell); 
		diagnostic->print_xs(xs, cell);
	    }

	// determine limiting event
	if (d_scatter < d_boundary && d_scatter < d_census)
	{
	    descriptor = "scatter";
	    dist_stream = d_scatter;
	}	
	else if (d_boundary < d_scatter && d_boundary < d_census)
	{
	    descriptor = "boundary";
	    dist_stream = d_boundary;
	}
	else 
	{
	    descriptor = "census";
	    dist_stream = d_census;
	}

	// IMC streaming
	stream_IMC(xs, tally, dist_stream);

	// scatter, effective or Thomson
	if (descriptor == "scatter")
	{
	    if (xs.get_sigma_thomson(cell) > 0.0)
	    {
		if (random.ran() < xs.get_sigeffscat(cell) /
		    (xs.get_sigeffscat(cell) + xs.get_sigma_thomson(cell)))
		    descriptor = "eff_scatter";
		else
		    descriptor = "thom_scatter";
	    }
	    else
		descriptor = "eff_scatter";

	    if (descriptor == "eff_scatter")
		tally.accum_n_effscat();
	    else if (descriptor == "thom_scatter")
		tally.accum_n_thomscat();

	    // accumulate momentum from before the scatter
	    tally.accumulate_momentum(cell, ew, omega);

	    // scatter the particle -- update direction cosines
	    scatter( mesh );

	    // accumulate momentum from after the scatter
	    tally.accumulate_momentum(cell, -ew, omega);
	}

	if (descriptor == "boundary")
        {
	    tally.accum_n_bndcross();
	    alive = surface(mesh, face);

	    if (descriptor == "reflection")
		tally.accum_n_reflections();

	    if (descriptor == "escape")
	    {
		tally.accum_n_escaped();
		tally.accum_ew_escaped(ew);
	    }

	}

	if (descriptor == "census")
	{
	    tally.accumulate_cen_info( cell, ew );
	    alive = false;
	}

	// do diagnostic print
	if (diagnostic)
	    diagnostic->print(*this);
    } 

    // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
// private transport member functions
//---------------------------------------------------------------------------//
// calculate everything about a collision

template<class MT>
bool Particle<MT>::collide(const MT &mesh, const Opacity<MT> &xs)
{   
    // status from collision
    bool status;

    // determine absorption or collision
    if (random.ran() <= xs.get_sigma_abs(cell) / xs.get_sigma_abs(cell))
    {
	descriptor = "absorption";
        status = false;
    }
    else
    {
        status = true;
        
	// calculate theta and phi (isotropic)
        double costheta, phi;
        costheta = 1 - 2 * random.ran();
        phi      = 2 * rtt_mc::global::pi * random.ran();

	// get new direction cosines
        mesh.get_Coord().calc_omega(costheta, phi, omega);
    }

    // return outcome of the event
    return status;
}

//---------------------------------------------------------------------------//
// perform an isotropic effective scatter 

template<class MT>
void Particle<MT>::scatter(const MT &mesh)
{   
    // calculate theta and phi (isotropic)
    double costheta, phi;
    costheta = 1 - 2 * random.ran();
    phi      = 2 * rtt_mc::global::pi * random.ran();
    
    // get new direction cosines
    mesh.get_Coord().calc_omega(costheta, phi, omega);
}


//---------------------------------------------------------------------------//
// do surface crossings

template<class MT>
bool Particle<MT>::surface(const MT &mesh, int face)
{
    // handle particles at a surface

    // status from surface crossing
    bool status;

    // determine the next cell
    int next_cell = mesh.next_cell(cell, face);

    // determine descriptor and outcome of this event

    if (next_cell == cell)
    {
	// reflection
	descriptor            = "reflection";
	vector<double> normal = mesh.get_normal(cell, face);
	double factor         = rtt_mc::global::dot(omega, normal);
	for (int i = 0; i < mesh.get_Coord().get_sdim(); i++)
	    omega[i] -= 2 * factor * normal[i];
	cell = next_cell;
    }
    else if (next_cell == 0)
    {
	// escape
	descriptor = "escape";
	cell       = next_cell;
    }
    else if (next_cell < 0)
    {
	// domain boundary crossing
	descriptor = "cross_boundary";
	cell       = next_cell;
    }
    else 
    {
	// continue streaming
	descriptor = "stream";
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
// diagnostic member functions
//---------------------------------------------------------------------------//
// print out a particle to some stream

template<class MT>
void Particle<MT>::print(ostream &output) const
{
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
// overloaded operators
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

    if (random.get_num() != rhs.random.get_num())
	return false;

    // if all these things check out then the particles are equal
    return true;
}

//===========================================================================//
// class Particle<MT>::Diagnostic
//===========================================================================//

//---------------------------------------------------------------------------//
// public diagnostic member functions
//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print(const Particle<MT> &particle)  const
{
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
    // print active particle (alive = true)
    output << " -- Particle is alive -- " << endl;

    // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << particle.descriptor.c_str() << endl;
    
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
    // print dead particle (alive = false)
    output << " -- Particle is dead -- " << endl;

    // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << particle.descriptor.c_str() << endl;
    
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
