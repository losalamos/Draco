//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Particle.t.hh"
 * \author Thomas M. Evans, Todd J. Urbatsch and Mike Buksas
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
	if (descriptor == VOL_EMISSION)
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
	    descriptor = SCATTER;
	    dist_stream = d_scatter;
	}	
	else if (d_boundary < d_scatter && d_boundary < d_census)
	{
	    descriptor = BOUNDARY;
	    dist_stream = d_boundary;
	}
	else 
	{
	    descriptor = CENSUS;
	    dist_stream = d_census;
	}

	// IMC streaming
	stream_IMC(xs, tally, dist_stream);

	// scatter, effective or Thomson
	if (descriptor == SCATTER)
	{
	    if (xs.get_sigma_thomson(cell) > 0.0)
	    {
		if (random.ran() < xs.get_sigeffscat(cell) /
		    (xs.get_sigeffscat(cell) + xs.get_sigma_thomson(cell)))
		    descriptor = EFF_SCATTER;
		else
		    descriptor = THOM_SCATTER;
	    }
	    else
		descriptor = EFF_SCATTER;

	    if (descriptor == EFF_SCATTER)
		tally.accum_n_effscat();
	    else if (descriptor == THOM_SCATTER)
		tally.accum_n_thomscat();

	    // accumulate momentum from before the scatter
	    tally.accumulate_momentum(cell, ew, omega);

	    // scatter the particle -- update direction cosines
	    scatter( mesh );

	    // accumulate momentum from after the scatter
	    tally.accumulate_momentum(cell, -ew, omega);
	}

	if (descriptor == BOUNDARY)
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

	if (descriptor == CENSUS)
	{
	    tally.accumulate_cen_info( cell, ew );
	    alive = false;
	    Check(rtt_mc::global::soft_equiv(time_left, 0.0));
	    time_left = 0.0;
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
	descriptor = ABSORPTION;
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
	   << setw(12) << get_descriptor(particle.descriptor) << endl;
    
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
	   << setw(12) << get_descriptor(particle.descriptor) << endl;
    
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

template<typename MT>
Particle<MT>::SP_Pack Particle<MT>::pack() const {

    // Require conditions

    // Create a smart pointer to a new Particle::Pack 
    return Particle<MT>::SP_Pack(new Particle::Pack(*this));

}

//===========================================================================//
// PARTICLE::PACK CLASS DEFINITIONS
//===========================================================================//

template <typename MT> int Particle<MT>::Pack::double_data_size = 0;
template <typename MT> int Particle<MT>::Pack::int_data_size    = 0;
template <typename MT> int Particle<MT>::Pack::char_data_size   = 0;
template <typename MT> bool Particle<MT>::Pack::setup_done      = false;

template <typename MT>
void Particle<MT>::Pack::set_sizes(int dim, int rnd_size)
{
    Require(dim>0 && dim<4); 
    Require(rnd_size);

    double_data_size = dim+6;  // This is all of the data, incl time_left
    int_data_size    = 2;
    char_data_size   = rnd_size;
    setup_done       = true;
}
    

/*!
 * \brief Setup function

 * Setup the sizes of the double, int and char data based on information from
 * the mesh and Random_Control objects. Should be called before constructing
 * Particle::Pack objects.  These methods are called by the Particle_Buffer
 * constructors

 */

template <typename MT> 
void Particle<MT>::Pack::setup_buffer_sizes(const MT &mesh, 
					    const rtt_rng::Rnd_Control &rcon)
{
    int dim = mesh.get_Coord().get_dim();  
    Check(dim>0 && dim<4);

    int rnd = rcon.get_size();
    Check(char_data_size <= rtt_rng::max_buffer);

    set_sizes(dim,rnd);
}


template <typename MT> 
void Particle<MT>::Pack::setup_buffer_sizes(const int dim,
					    const rtt_rng::Rnd_Control &rcon)
{
    Require(dim>0 && dim<4);

    int rnd = rcon.get_size();
    Check(char_data_size <= rtt_rng::max_buffer);

    set_sizes(dim,rnd);
}


/*!
 * \brief Constructor.

 * Construct a Particle::Pack instance. Once allocated char, double and int
 * data is given to the constructor via pointers, the Pack object owns it and 
 * is resonsible for deallocation. Normally, Pack objects are created via
 * Particle::pack() which returns a smart pointer to the Pack object.

 */

template <typename MT>
Particle<MT>::Pack::Pack(double* d, int* i, char* c) :
    double_data(d), int_data(i), char_data(c) 
{ /* Assume control of the given pointers */ }


/*!
 * \brief Constructor.
 * Construct a Particle::Pack instance from a particle.
 */
template <typename MT>
Particle<MT>::Pack::Pack(const Particle<MT>& particle)
{
    Require(particle.cell);

    // Temp storage
    char* bytes;
    
    // Set data sizes
    int dim = particle.r.size();    Check (dim>0 && dim<4);

    // Compute sizes locally and compare. Note that this duplicates the size info
    int d_size = dim + 6; 
    int i_size = 2;
    int c_size = pack_sprng(particle.random.get_id(), &bytes); Check(bytes);

    // if this is the first construction, use these values to set up the
    // static data, else, check for consistency
    if (!setup_done) set_sizes(dim, c_size);
    else {
	Check(d_size == double_data_size);
	Check(i_size == int_data_size);
	Check(c_size == char_data_size);
    }

    // Allocate the storage
    double_data = new double[double_data_size];
    int_data    = new int   [int_data_size];
    char_data   = new char  [char_data_size]; 

    // Insert the values
    // double data
    double_data[0] = particle.time_left;
    double_data[1] = particle.ew;
    double_data[2] = particle.fraction;
    std::copy(particle.omega.begin(), particle.omega.end(), double_data+3);
    std::copy(particle.r.begin(),     particle.r.end(),     double_data+6);

    // int data
    int_data[0] = particle.cell;
    int_data[1] = particle.random.get_num();

    // char data
    std::copy(bytes, bytes+char_data_size, char_data);
    std::free(bytes);

}


//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.

 * Do copy construction while preserving memory.  This is not a reference 
 * counted class so data is copied from one class to the other during
 * function calls and the like (wherever a copy constructor is called).

 */

template <typename MT>
Particle<MT>::Pack::Pack(const Pack &rhs) : 
    double_data(new double[rhs.double_data_size]), 
    int_data(new int[rhs.int_data_size]),
    char_data(new char[rhs.char_data_size])
{
    // copy all three data arrays
    std::copy(rhs.double_data, rhs.double_data+double_data_size,   double_data);
    std::copy(rhs.int_data,    rhs.int_data   +int_data_size,      int_data);
    std::copy(rhs.char_data,   rhs.char_data  +rhs.char_data_size, char_data);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.

 * Cleans up memory when the Pack object goes out of scope.  Once allocated
 * pointers are given to the Pack object the Pack object takes control of
 * them.

 */

template <typename MT>
Particle<MT>::Pack::~Pack()
{
    delete [] double_data;
    delete [] int_data;
    delete [] char_data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack the Particle

 * Unpacks and returns a smart pointer to a new Particle.

 * \return smart pointer to the unpacked particle

 */

template <typename MT>
Particle<MT>::SP_Particle Particle<MT>::Pack::unpack() const
{
    // calls the static version of unpack with it's own data
    return unpack(double_data_size, double_data,
		  int_data_size,    int_data,
		  char_data_size,   char_data);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Unpack the Particle represented as raw data

 * Unpacks and returns a smart pointer to a new Particle.

 * \return smart pointer to the unpacked particle

 */

template <typename MT>
Particle<MT>::SP_Particle Particle<MT>::Pack::unpack(
    int double_size, double* double_data, 
    int int_size,    int*    int_data,
    int char_size,   char*   char_data)
{
    Require(double_size>6);  Require(double_data);
    Require(int_size==2);    Require(int_data);
    Require(char_size>0);    Require(char_data);

    // compute the dimension of the space based on the size of the data:
    int dim = double_size - 6;
    Check(dim>0 && dim<4);

    // double data
    double t_left = double_data[0];  Check(t_left>=0.0);
    double ew     = double_data[1];  Check(ew>0.0);
    double frac   = double_data[2];  Check(frac>0.0);
    vector<double> omega(double_data+3, double_data+6);
    vector<double> r    (double_data+6, double_data+6+dim);

    // int data
    int cell = int_data[0];  Check(cell);
    int rnum = int_data[1];  Check(rnum>=0);

    // char data
    char *bytes = new char[char_size];
    std::copy(char_data, char_data+char_size, bytes);

    // random object, from int and char data
    int *rnid = unpack_sprng(bytes);
    rtt_rng::Sprng random(rnid, rnum);
    delete [] bytes;
    
    return SP_Particle (new Particle(r, omega, ew, cell, random, frac, t_left));

}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Particle.t.hh
//---------------------------------------------------------------------------//
