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

#ifndef __imc_Particle_hh__
#define __imc_Particle_hh__

#include "Opacity.hh"
#include "Tally.hh"
#include "Global.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

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
 * rtt_imc::Particle_Buffer.

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
// 
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
		     SCATTER=100, 
		     LOW_WEIGHT=101,
		     EFF_SCATTER=102, 
		     THOM_SCATTER=103, 
		     REFLECTION=200, 
		     STREAM=201, 
		     ESCAPE=202, 
		     CROSS_BOUNDARY=203, 
		     ABSORPTION=204,
		     BOUNDARY=205, 
		     CENSUS=300, 
		     KILLED=1000};

    // Forward declaration of Pack class
    struct Pack;

    // Public typedefs:
    typedef rtt_dsxx::SP<Particle::Pack> SP_Pack;
    typedef rtt_dsxx::SP<Particle>       SP_Particle;
    typedef std::vector<double>          sf_double;
    typedef rtt_rng::Sprng               Rnd_Type;
    typedef std::string                  std_string;

    /*!
     * \class Particle::Diagnostic
     * \brief Diagnostic class for tracking particle histories.
     */
    class Diagnostic
    {
      private:
	// Stream output is sent to.
	std::ostream &output;

	// Boolean for detailed diagnostic.
	bool detail;

      public:
	//! Constructor.
	Diagnostic(std::ostream &output_, bool detail_ = false) 
	    : output(output_), detail(detail_) {}

	//! Access detail switch (true gives more destail in printouts). 
	bool detail_status() const { return detail; }

	// Diagnostic print functions.
	void print(const Particle<MT> &) const;
	void print_alive(const Particle<MT> &) const;
	void print_dead(const Particle<MT> &) const;
	void print_dist(double, double, double, int) const;
	void print_xs(const Opacity<MT> &, int) const;

	//! Inline output formatters
	inline void header() const;
    };

    // Friend declarations.
    friend class Diagnostic;
    friend class Pack;

  private:
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
    Rnd_Type random;

  private:
    // >>> IMPLEMENTATION

    // Stream a distance d.
    inline void stream(double);  

    // Stream a distance d.
    inline void stream_IMC(const Opacity<MT> &, Tally<MT> &, double);

    // Collision, return a false if particle is absorbed.
    bool collide(const MT &, const Opacity<MT> &);

    // Effective scatter.
    void scatter(const MT & );

    // Surface crossings, return a false if particle escapes.
    bool surface(const MT &, int);

    // IMC transport step.
    void trans_IMC(const MT &, const Opacity<MT> &);
    
    // DDMC transport step.
    void trans_DDMC(const MT &, const Opacity<MT> &);

  public:
    // Particle constructor.
    inline Particle(const sf_double &, const sf_double &, double, int,
		    Rnd_Type, double = 1, double = 1, int = BORN);

    // >>> TRANSPORT INTERFACE

    // IMC transport step.
    void transport(const MT &, const Opacity<MT> &, Tally<MT> &,
		   rtt_dsxx::SP<Diagnostic> = rtt_dsxx::SP<Diagnostic>()); 

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
    const Rnd_Type& get_random() const { return random; }

    //! Get the particle descriptor.
    int get_descriptor() const { return descriptor; }

    //! Convert an integer descriptor into a string.
    inline static std_string get_descriptor(int);
    
    //! Convert a string descriptor into an integer.
    inline static int get_index(std_string);

    // >>> SET FUNCTIONS
    
    //! Reset the particle status to alive (true).
    void reset_status() { alive = true; }
    
    //! Set the particle status to dead (false).
    void kill_particle() { alive = false; }

    //! Set the particle's random number.
    void set_random(const Rnd_Type &ran) { random = ran; }

    //! Set the particle's time-left fraction.
    void set_time_left(double t) { time_left = t; }

    //! Set the particle's descriptor.
    void set_descriptor(int desc) { descriptor = desc; }

    //! Set the particle's energy weight.
    void set_ew(double new_ew) { ew = new_ew; }

    //! Set the particle's cell index.
    void set_cell(int new_cell) { cell = new_cell; }

    // >>> DIAGNOSTIC FUNCTIONS

    // Public diagnostic services.
    void print(std::ostream &) const;

    // Overloaded equals operator.
    bool operator==(const Particle<MT> &) const;

    //! Overloaded not equals operator.
    bool operator!=(const Particle<MT> &p) const { return !(*this == p); } 

    // Pack function
    SP_Pack pack() const;

};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// output ascii version of Particle

template<class MT>
inline std::ostream& operator<<(std::ostream &output, 
				const Particle<MT> &object)
{
    object.print(output);
    return output;
}

//---------------------------------------------------------------------------//
// inline functions for Particle<MT>::Diagnostic
//---------------------------------------------------------------------------//
// Particle<MT>::Diagnostic inline functions

template<class MT>
inline void Particle<MT>::Diagnostic::header() const 
{ 
    output << "*** PARTICLE HISTORY ***" << std::endl; 
    output << "------------------------" << std::endl;
}

//---------------------------------------------------------------------------//
// Particle<MT> inline functions
//---------------------------------------------------------------------------//
// Particle<MT> constructor

template<class MT>
inline Particle<MT>::Particle(const sf_double& r_, 
			      const sf_double& omega_, 
			      double ew_,
			      int cell_, 
			      Rnd_Type random_, 
			      double frac,
			      double tleft, 
			      int desc)
    : ew(ew_), 
      r(r_),
      omega(omega_),
      cell(cell_),
      time_left(tleft), 
      fraction(frac), 
      alive(true), 
      descriptor(desc),
      random(random_)
{
    // non-default particle constructor
    Ensure (r.size() < 4);
    Ensure (omega.size() == 3);
    Ensure (cell > 0);
}

//---------------------------------------------------------------------------//

template<class MT>
inline void Particle<MT>::stream(double distance)
{
    // calculate new location when Particle streams
    for (int i = 0; i <= r.size()-1; i++)
	r[i] = r[i] + distance * omega[i];
}

//---------------------------------------------------------------------------//

template<class MT>
inline void Particle<MT>::stream_IMC(const Opacity<MT> &xs, Tally<MT> &tally,
				     double distance)
{
    // hardwire minimum energy weight fraction; limit distance accordingly
    double minwt_frac = 0.01;
    double min_arg    = log(0.1 * minwt_frac);
    double argument = -xs.get_sigeffabs(cell) * distance;
    if (argument < min_arg) argument = min_arg;

    // calculate multiplicative reduction in energy-weight; 
    // calculate new energy weight; change in energy-weight
    double factor = exp(argument);
    double new_ew = ew * factor;
    double del_ew = ew - new_ew;

    // accumulate tallies from time-rate absorption.  ewpl ==
    // energy-weighted-path-length = int_0^d e^(-sig*x) dx =
    // (1/sig)ew(1-e^(-sig*x)), or if sig=0, ewpl = ew*d.
    tally.deposit_energy( cell, del_ew );
    tally.accumulate_momentum(cell, del_ew, omega);
    if (xs.get_sigeffabs(cell) > 0)
	tally.accumulate_ewpl(cell, del_ew / xs.get_sigeffabs(cell) );
    else if (xs.get_sigeffabs(cell) == 0)
	tally.accumulate_ewpl(cell, distance * ew);
    else if (xs.get_sigeffabs(cell) < 0)
	Insist (0, "Effective absorption is negative!");

    // update the fraction of the particle's original weight
    fraction *= factor;

    // kill particle and deposit its energy if its fraction is too low;
    // or update particle energy-weight, time_left, and position.
    if (fraction < minwt_frac)
    {
	tally.deposit_energy( cell, new_ew );
	tally.accum_n_killed();
	tally.accum_ew_killed( new_ew );
	tally.accumulate_momentum(cell, new_ew, omega);
	descriptor = KILLED;
	alive = false;
    }
    else
    {
	ew = new_ew;
	time_left -= distance / rtt_mc::global::c;
	for (int i = 0; i <= r.size()-1; i++)
	    r[i] = r[i] + distance * omega[i];
    }
}

//---------------------------------------------------------------------------//
// convert a particle event descriptor string into an int

template<class MT>
inline int Particle<MT>::get_index(std::string desc)
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

/* Changed */

template<class MT>
inline std::string Particle<MT>::get_descriptor(int index)
{

    switch (index) {

    // born descriptors
    case BORN: 	         return "born";
    case CENSUS_BORN:	 return "census_born";
    case BOUNDARY_BORN:	 return "boundary_born";
    case VOL_EMISSION:	 return "vol_emission";
    case SURFACE_SOURCE: return "surface_source";

    // collision event descriptors
    case SCATTER:        return "scatter";
    case LOW_WEIGHT: 	 return "low_weight";
    case EFF_SCATTER:	 return "eff_scatter";
    case THOM_SCATTER:	 return "thom_scatter";

    // streaming descriptors
    case REFLECTION:	 return "reflection";
    case STREAM:	 return "stream";
    case ESCAPE:	 return "escape";
    case CROSS_BOUNDARY: return "cross_boundary";

    // time and census descriptors
    case CENSUS:	 return "census";

    // death
    case KILLED:	 return "killed";

    default:
	Insist(0,"Unrecognized descriptor number");
	return "";
	    
    }    

}


//===========================================================================//
/*!  
 * \struct Particle::Pack 
 
 * \brief Nested class for packing particles into raw data for writing
 * or communication.
 
 */
//===========================================================================//

template <typename MT>
struct Particle<MT>::Pack { 
    
  private:
    
    // Size info is static for class wide consistency and access
    static int double_data_size;
    static int int_data_size;
    static int char_data_size;
    static bool setup_done;

    double *double_data; 
    int *int_data;       
    char*char_data;      

    // Disallow assignment
    const Pack& operator=(const Pack &);

    // Internal size set routine
    static void set_sizes(int, int);

  public:
    
    /* Setup */

    static void setup_buffer_sizes(const MT& mesh, 
				   const rtt_rng::Rnd_Control &rcon);

    static void setup_buffer_sizes(const int dim,
				   const rtt_rng::Rnd_Control &rcon);

    /* Structors */
    
    // Construct from three pointers
    Pack(double*, int*, char*); 

    // Construct from Particle
    Pack(const Particle&);
    
    // Copy
    Pack(const Pack&);
    
    // Destroy
    ~Pack();
    
    /* Accessors */
    
    //! Get const pointer to beginning of integer data stream
    const int* int_begin() const {return int_data; }
    
    //! Get const pointer to begining of char data stream
    const char* char_begin() const {return char_data; }
    
    //! Get const pointer to begining of double double stream
    const double* double_begin() const {return double_data; }
    
    //! Get const pointer to one past end of integer data stream
    const int* int_end() const {return int_data+int_data_size; }
    
    //! Get const pointer to one past end of char data stream
    const char* char_end() const {return char_data+char_data_size; }
    
    //! Get const pointer to one past end of double data stream
    const double* double_end() const {return double_data+double_data_size; }
    
    //! Unpack function
    SP_Particle unpack() const;

    /* Functions for size information is static to assist in constructing
       Particle_Buffer objects */

    //! Get size of integer data stream
    static int get_int_size() { return int_data_size; }
    
    //! Get size of char data stream
    static int get_char_size() { return char_data_size; }
    
    //! Get size of double data stream
    static int get_double_size() { return double_data_size; }

    //! Get setup status of class
    static bool get_setup() { return setup_done; }
    
    //! Unpack from provided data
    static SP_Particle unpack(int, double*, int, int*, int, char*);
    
}; // End of class Particle::Pack


} // end namespace rtt_imc

#endif                          // __imc_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Particle.hh
//---------------------------------------------------------------------------//
