//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Particle.hh 
 * \author Thomas M. Evans, Todd J. Urbatsch and Mike Buksas
 * \date   Fri Jan 30 17:04:24 1998
 * \brief  Particle class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Particle_hh
#define rtt_imc_Particle_hh

#include "Opacity.hh"
#include "Tally.hh"
#include "Global.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"
#include <vector>
#include <string>
#include <iostream> 
#include <cmath>
#include <iomanip>
#include <algorithm>

namespace rtt_imc 
{

//===========================================================================//
/*!
 * \class Particle

 * \brief Defines IMC Fleck and Cummings particles.

 * The particle class defines the attributes, behavior, and relationships for
 * Fleck and Cummings IMC particles.  The primary particle attributes are
 * position, direction, weight, random number object, cell index, time, and
 * fraction.

 * The particles defined by this class carry along their own random number
 * object.  Thus, reproducible Monte Carlo calculations may be performed.

 * Particle containers and buffers (for communication) are defined in
 * rtt_mc::Particle_Buffer and rtt_mc::Particle_Stack.  Input/Output
 * functions are in rtt_mc::Particle_IO.  A communication class for parallel
 * transport (that uses the Particle_Buffer classes) is defined in the
 * rtt_mc::Communicator class.

 * The primary operation for IMC particles is transport().  Given a mesh,
 * opacity, and material data, the particles are capable of transporting
 * until some termination point.  The standard termination points are:

 * \arg \b escape exiting a problem boundary

 * \arg \b cross_boundary crossing a processor boundary

 * \arg \b census going to census

 * \arg \b killed reaching a weight cutoff

 * If any of these termination points are reached the particle ceases
 * transport.  The rtt_imc::Tally class tallies all of the event data that
 * the particle undergoes during a transport step.

 * The class also includes several accessors and modifiers for adapting a
 * particle before or after transport.  The status of the particle can be
 * accessed through the status() function.  Particles must have true status
 * before they can be transported.  The reset_status() function allows a user
 * to "turn on" a particle.

*/
/*!
 * \example imc/test/tstParticle.cc

 * Example usage of the rtt_imc::Particle and rtt_imc::Particle_Buffer
 * classes.

*/
// revision history:
// -----------------
//  0) original
//  1)   2-3-98 : made Particle a class template parameterized on Mesh_type
//  2)   4-3-98 : added Particle-Stack struct; added friendship to list and 
//                stack stl template classes.  If add stacks based on other
//                stl containers don't forget to add friendship in Particle 
//  3)   4-6-98 : moved the default constructor to public with an assert(0); 
//                this way, if anyone tries to use it the code will crash at 
//                runtime, yet, the STL containers still can see the 
//                constructor to make the containers, can't give objects to 
//                initialize the containers because there are no conversion 
//                constructors in Particle; made Particle_Stack parameterized 
//                on PT (Particle-type)
//  4)   5-5-98 : removed source() function and parameterization on random
//                number type, updated contructor to work with new source
//                classes, added dump and changed random number to Sprng 
//                types.
//  5)  5-12-98 : added Particle_Buffer nested class to aid in persistence of 
//                particles
//  6)  5-26-98 : added static functions for getting the particle descriptor
//                and converting it to an int
//  7)   7-6-98 : added functions to return particle information to needed
//                services 
//  8)   7-9-98 : added set functions (ew and random) for combing
//  9)  7-28-98 : added Thomson scattering
// 10)  9-20-98 : added overloaded operator== to check equality of Particle
//                data for testing, this only checks if the Particle's random
//                number has the SAME streamnum, not the same ID, because you
//                will get the same stream for the same streamnum even if you
//                have a different pointer pointing to the state data
// 11)  4-29-99 : removed using declarations from namespace
// 12) 18-AUG-00: updated surface() mesh.next_cell() to include the position
//                vector for AMR support. 
// 13) 28-AUG-00: modified the energy-weighted path length to account for a
//                zero opacity.  
// 14) 27-APR-01: added pack struct and member function
// 15) 30-MAY-01: added Descriptor enum to represent particle state, removed
//                string descriptor. See memo CCS-4:01-19(U)
// 16) 26-JUL-01: cleaned up file format; made get_descriptor return the
//                enumeration value (integer)
// 17) 24-SEP-01: switched to hybrid aborption behavior and added streaming
//                to cutoff. See memo CCS-4:01-33(U)
// 18) 20-DEC-01: redesigned the pack infrastructure
// 19) 09-JAN-02: updated particle to work with new Particle_Buffer; added
//                packed size calculation static function
// 20) 29-JAN-02: added Gray_Particle and Multigroup_Particle
//===========================================================================//

template<class MT>
class Particle
{
  public: 
    // >>> NESTED TYPES

    // Descriptor enumeration for particle events.
    enum Descriptor {BORN=1, 
		     CENSUS_BORN=2, 
		     BOUNDARY_BORN=3, 
		     VOL_EMISSION=4,
		     SURFACE_SOURCE=5,
		     UNPACKED=6,
		     SCATTER=100, 
		     LOW_WEIGHT=101,
		     EFF_SCATTER=102, 
		     THOM_SCATTER=103, 
		     COLLISION=104, 
		     CUTOFF=105,
		     REFLECTION=200, 
		     STREAM=201, 
		     ESCAPE=202, 
		     CROSS_BOUNDARY=203, 
		     ABSORPTION=204,
		     BOUNDARY=205, 
		     CENSUS=300, 
		     RANDOM_WALK=400,
		     DDIMC=401,
		     KILLED=1000};

    /*!
     * \brief Base Diagnostic class for tracking particle histories.
     */
    class Diagnostic
    {
      protected:
	// Stream output is sent to.
	std::ostream &output;
    
	// Boolean for detailed diagnostic.
	bool detail;
    
      public:
	//! Constructor.
	Diagnostic(std::ostream &output_, bool detail_) 
	    : output(output_), detail(detail_) {}

	//! Destructor for inheritance chain.
	virtual ~Diagnostic() {/*...*/}
    
	//! Access detail switch (true gives more destail in printouts). 
	bool detail_status() const { return detail; }
    
	// Diagnostic print functions.
	virtual void print(const Particle<MT> &) const;
	virtual void print_alive(const Particle<MT> &) const;
	virtual void print_dead(const Particle<MT> &) const;
	virtual void print_dist(double, double, double, int) const;
    
	//! Inline output formatters
	inline void header() const;
    };

    // Useful typedefs.
    typedef std::vector<double>      sf_double;
    typedef rtt_rng::Sprng           Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>   SP_Rnd_Type;
    typedef std::string              std_string;
    
    // friend declarations
    friend class Diagnostic;

  protected:

    // >> PARAMERERS:

    // Energy-weight cutoff between analog and implicit absorption
    static const double minwt_frac;

    // >>> DATA

    // Particle energy-weight.
    double ew;

    // Particle location.
    sf_double r;

    // Particle direction.
    sf_double omega;

    // Particle cell.
    int cell;

    // Time remaining in this time step.
    double time_left;
    
    // Fraction of original energy weight.
    double fraction;

    // Status of particle.
    bool alive;
    
    // Event type descriptor.
    int descriptor;

    // Random number object.
    SP_Rnd_Type random;

  protected:

    // >>> IMPLEMENTATION

    // Stream a distance d.
    inline void stream(double);  

    // Perform streaming operations specific to particles undergoing analog
    // absorption
    inline void stream_analog_capture  (Tally<MT> &, double);

    // Surface crossings, return a false if particle escapes.
    inline bool surface(const MT &, int);

    //! Check for energy-weight fraction below Implicit/Analog cutoff. 
    bool use_analog_absorption() { return (fraction <= minwt_frac); }

    // Process a particle going into the census
    inline void census_event(Tally<MT> &);

    // Process a particle crossing a boundary
    inline void boundary_event(const MT&, Tally<MT>&, int);

    // Perform an isotropic scatter.
    inline void scatter(const MT &);

  public:
    // Particle constructor.
    inline Particle(const sf_double &, const sf_double &, double, int,
		    Rnd_Type, double = 1, double = 1, int = BORN);

    //! Default constructor.
    Particle() {/*...*/}

    //! Virtual destructor.
    virtual inline ~Particle() = 0;
		      
    // >>> ACCESSORS

    //! Return the particle status.
    bool status() const { return alive; }

    //! Get the particle's current cell index.
    int get_cell() const { return cell; }

    //! Get the particle's energy weight.
    double get_ew() const { return ew; }

    //! Get the particle's position.
    const sf_double& get_r() const { return r; }

    //! Get the particle's direction.
    const sf_double& get_omega() const { return omega; }

    //! Get the particle's random number object.
    const Rnd_Type& get_random() const { return *random; }

    //! Get the particle descriptor.
    int get_descriptor() const { return descriptor; }

    //! Convert an integer descriptor into a string.
    static std_string convert_descriptor_to_string(int);
    
    //! Convert a string descriptor into an integer.
    static int convert_string_to_descriptor(std_string);

    // >>> SET FUNCTIONS
    
    //! Reset the particle status to alive (true).
    void reset_status() { alive = true; }
    
    //! Set the particle status to dead (false).
    void kill_particle() { alive = false; }

    //! Set the particle's random number.
    void set_random(const Rnd_Type &ran) { random = new Rnd_Type(ran); }

    //! Set the particle's time-left fraction.
    void set_time_left(double t) { time_left = t; }

    //! Set the particle's descriptor.
    void set_descriptor(int desc) { descriptor = desc; }

    //! Set the particle's energy weight.
    void set_ew(double new_ew) { ew = new_ew; }

    //! Set the particle's cell index.
    void set_cell(int new_cell) { cell = new_cell; }

    // >>> DIAGNOSTIC FUNCTIONS
    virtual void print(std::ostream &) const;
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Output a particle.
 */
template<class MT>
std::ostream& operator<<(std::ostream       &output, 
			 const Particle<MT> &object)
{
    object.print(output);
    return output;
}

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Base class Particle constructor.
 *
 * The constructor is declared inline for optimization purposes.
 */
template<class MT>
Particle<MT>::Particle(const sf_double &r_, 
		       const sf_double &omega_, 
		       double           ew_,
		       int              cell_, 
		       Rnd_Type         random_, 
		       double           frac,
		       double           tleft, 
		       int              desc)
    : ew(ew_), 
      r(r_),
      omega(omega_),
      cell(cell_),
      time_left(tleft), 
      fraction(frac), 
      alive(true), 
      descriptor(desc),
      random(new Rnd_Type(random_))
{
    // non-default particle constructor
    Check (r.size() < 4);
    Check (omega.size() == 3);
    Check (cell > 0);
}

//---------------------------------------------------------------------------//

template<class MT>
Particle<MT>::~Particle()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Move particle a distance.
 */
template<class MT>
void Particle<MT>::stream(double distance)
{
    // calculate new location when Particle streams
    for (int i = 0; i <= r.size()-1; i++)
	r[i] = r[i] + distance * omega[i];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream the particle during analog transport.
 */
template<class MT>
void Particle<MT>::stream_analog_capture(Tally<MT> &tally, double distance)
{
    Check(distance >= 0);

    // accumulate the energy-weighted pathlength
    tally.accumulate_ewpl(cell, distance * ew);

    // Physically transport the particle
    stream(distance); 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle going into census.
 */
template<class MT>
void Particle<MT>::census_event(Tally<MT> &tally)
{
    // tally census data about the particle
    tally.accumulate_cen_info( cell, ew );
    
    // set status to dead
    alive = false;

    Check(rtt_mc::global::soft_equiv(time_left, 0.0));

    // set time left to identically zero
    time_left = 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a boundary crossing event.
 */
template<class MT>
void Particle<MT>::boundary_event(const MT  &mesh, 
				  Tally<MT> &tally, 
				  int        face)
{
    // tally the boundary crossing
    tally.accum_n_bndcross();

    // determine the status at the surface
    alive = surface(mesh, face);
    
    // tally the reflection
    if (descriptor == REFLECTION)
	tally.accum_n_reflections();
    
    // tally the escape
    if (descriptor == ESCAPE)
    {
	tally.accum_n_escaped();
	tally.accum_ew_escaped(ew);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Handle surface crossing.
 *
 * \return false if particle leaves system
 */
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
	Check (normal.size() == 3);
	
	// specularly reflect angle
	double factor = rtt_mc::global::dot(omega, normal);
	for (int i = 0; i < 3; i++)
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
/*!
 * \brief Perform an isotropic scatter.
 */
template<class MT>
void Particle<MT>::scatter(const MT &mesh)
{   
    // calculate theta and phi (isotropic)
    double costheta = 1.0 - 2.0 * random->ran();
    double phi      = 2.0 * rtt_mc::global::pi * random->ran();
    
    // get new direction cosines
    mesh.get_Coord().calc_omega(costheta, phi, omega);
}

//---------------------------------------------------------------------------//
// INCLUDE Particle.t.hh SO AUTOMATIC INSTANTIATION WILL WORK 
//---------------------------------------------------------------------------//

#include "Particle.t.hh"

} // end namespace rtt_imc

#endif                          // rtt_imc_Particle_hh

//---------------------------------------------------------------------------//
//                              end of imc/Particle.hh
//---------------------------------------------------------------------------//
