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
#include "ds++/Packing_Utils.hh"
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
		     KILLED=1000};

    // Useful typedefs.
    typedef std::vector<double>      sf_double;
    typedef rtt_rng::Sprng           Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>   SP_Rnd_Type;
    typedef std::string              std_string;

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

  private:

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

  private:
    // >>> IMPLEMENTATION

    // Stream a distance d.
    inline void stream(double);  

    //! Perform streaming operations specific to particles undergoing
    //implicit absorption
    inline void stream_implicit_capture(const Opacity<MT> &, Tally<MT> &, 
					double);

    //! Perform streaming operations specific to particles undergoing analog
    // absorption
    inline void stream_analog_capture  (Tally<MT> &, double);

    // Effective scatter.
    void scatter(const MT & );

    // Surface crossings, return a false if particle escapes.
    bool surface(const MT &, int);

    //! Check for energy-weight fraction below Implicit/Analog cutoff. 
    bool use_analog_absorption() { return (fraction <= minwt_frac); }

    //! Process a combined scattering/collision event
    void collision_event(const MT&, Tally<MT>&, double, double, double);

    //! Process a particle going into the census
    void census_event(Tally<MT> &);

    //! Process a particle crossing a boundary
    void boundary_event(const MT&, Tally<MT>&, int);

  public:
    // Particle constructor.
    inline Particle(const sf_double &, const sf_double &, double, int,
		    Rnd_Type, double = 1, double = 1, int = BORN);

    // Unpacking constructor.
    inline Particle(const std::vector<char> &);

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

    // Public diagnostic services.
    void print(std::ostream &) const;

    // Overloaded equals operator.
    bool operator==(const Particle<MT> &) const;

    //! Overloaded not equals operator.
    bool operator!=(const Particle<MT> &p) const { return !(*this == p); } 

    // Pack function
    inline std::vector<char> pack() const;
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
// output ascii version of Particle

template<class MT>
std::ostream& operator<<(std::ostream &output, 
			 const Particle<MT> &object)
{
    object.print(output);
    return output;
}

//---------------------------------------------------------------------------//
// PARTICLE<MT> INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Regular particle constructor.
 *
 * The constructor is declared inline for optimization purposes.
 */
template<class MT>
Particle<MT>::Particle(const sf_double& r_, 
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
      random(new Rnd_Type(random_))
{
    // non-default particle constructor
    Ensure (r.size() < 4);
    Ensure (omega.size() == 3);
    Ensure (cell > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor is used to unpack a particle that has been packed with
 * the Particle::pack() function.
 *
 * It is declared inline for optimization.
 */
template<class MT>
Particle<MT>::Particle(const std::vector<char> &packed)
{
    Require (packed.size() >= 3 * sizeof(int) + 6 * sizeof(double));

    // make an unpacker
    rtt_dsxx::Unpacker u;
    
    // set it
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the spatial dimension of the particle
    int dimension = 0;
    u >> dimension;
    Check (dimension > 0 && dimension <= 3);

    // size the dimension and direction 
    r.resize(dimension);
    omega.resize(3);

    // unpack the position
    for (int i = 0; i < dimension; i++)
	u >> r[i];

    // unpack the rest of the data
    u >> omega[0] >> omega[1] >> omega[2] >> cell >> ew >> time_left
      >> fraction;
    Check (time_left >= 0.0);
    Check (fraction  >= 0.0);
    Check (cell      >  0);
    Check (ew        >= 0.0);

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
    random = new Rnd_Type(prn);
    Check (random->get_num() >= 0);
    Check (random->get_id());

    // assign the descriptor and status
    descriptor = UNPACKED;
    alive      = true;
    
    Ensure (status());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a particle into a char stream for communication and
 * persistence. 
 */
template<class MT>
std::vector<char> Particle<MT>::pack() const
{
    Require (omega.size() == 3);

    // make a packer
    rtt_dsxx::Packer p;

    // first pack the random number state
    std::vector<char> prn = random->pack();

    // determine the size of the packed particle: 1 int for cell, + 1 int for
    // size of packed RN state + 1 int for dimension of space; dimension +
    // 6 doubles; + size of RN state chars
    int size = 3 * sizeof(int) + (r.size() + 6) * sizeof(double) + prn.size();

    // set the packed buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the spatial dimension
    p << static_cast<int>(r.size());
    
    // pack the dimension
    for (int i = 0; i < r.size(); i++)
	p << r[i];
    
    // pack the rest of the data
    p << omega[0] << omega[1] << omega[2] << cell << ew << time_left
      << fraction;

    // pack the RN state
    p << static_cast<int>(prn.size());
    for (int i = 0; i < prn.size(); i++)
	p << prn[i];

    Ensure (p.get_ptr() == &packed[0] + size);
    return packed;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::stream(double distance)
{
    // calculate new location when Particle streams
    for (int i = 0; i <= r.size()-1; i++)
	r[i] = r[i] + distance * omega[i];
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::stream_implicit_capture(const Opacity<MT> &xs, 
					   Tally<MT> &tally,
					   double distance)
{
    Check(distance>=0)

    double argument = -xs.get_sigeffabs(cell) * distance;

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

    // update particle energy-weight
    ew = new_ew;

    Check(ew)

    // Physically transport the particle
    stream(distance); 

}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::stream_analog_capture(Tally<MT> &tally, 
						double distance)
{
    Check(distance>=0)

    tally.accumulate_ewpl(cell, distance * ew);

    // Physically transport the particle
    stream(distance); 
}

} // end namespace rtt_imc

#endif                          // __imc_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Particle.hh
//---------------------------------------------------------------------------//
