//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Particle.hh 
 * \author Thomas M. Evans and Todd J. Urbatsch
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
// 
//===========================================================================//

template<class MT>
class Particle
{
  public: 

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
    template<class PT> friend class Particle_Buffer;

  private:
    //>>> DATA
    // Particle energy-weight.
    double ew;

    // Particle location.
    std::vector<double> r;

    // Particle direction.
    std::vector<double> omega;

    // Particle cell.
    int cell;

    // Time remaining in this time step.
    double time_left;
    
    // Fraction of original energy weight.
    double fraction;

    // Status of particle.
    bool alive;
    
    // Event type descriptor.
    std::string descriptor;

    // Random number object.
    rtt_rng::Sprng random;

  private:
    //>>> IMPLEMENTATION

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
    inline Particle(std::vector<double>, std::vector<double>, double, int,
		    rtt_rng::Sprng, double = 1, double = 1, 
		    std::string = "born");

    // Null constructor required as kluge for the STL containers which need a
    // default constructor, this calls an assert(0) so you can't use it.
    Particle() { Insist (0, "You tried to default construct a Particle!"); }

    //>>> TRANSPORT INTERFACE

    // IMC transport step.
    void transport(const MT &, const Opacity<MT> &, Tally<MT> &,
		   rtt_dsxx::SP<Diagnostic> = rtt_dsxx::SP<Diagnostic>()); 

    //>>> SERVICES

    //! Return the particle status.
    bool status() const { return alive; }
    
    //! Reset the particle status to alive (true).
    void reset_status() { alive = true; }
    
    //! Set the particle status to dead (false).
    void kill_particle() { alive = false; }

    // Accessors.
    int get_cell() const { return cell; }
    double get_ew() const { return ew; }
    std::vector<double>   get_omega() const { return omega; }
    const rtt_rng::Sprng& get_random() const { return random; }

    // Set functions for sourcing particles.
    void set_random(rtt_rng::Sprng &ran) { random = ran; }
    void set_time_left(double t) { time_left = t; }
    void set_descriptor(std::string s) { descriptor = s; }
    void set_ew(double new_ew) { ew = new_ew; }
    void set_cell(int new_cell) { cell = new_cell; }

    // Transport descriptors.
    std::string desc() const { return descriptor; }
    inline static int get_index(std::string);
    inline static std::string get_descriptor(int);

    // Public diagnostic services.
    void print(std::ostream &) const;

    // Overloaded operators.
    bool operator==(const Particle<MT> &) const;
    bool operator!=(const Particle<MT> &p) const { return !(*this == p); } 
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
inline Particle<MT>::Particle(std::vector<double> r_, 
			      std::vector<double> omega_, 
			      double ew_, int cell_, 
			      rtt_rng::Sprng random_, 
			      double frac, double tleft, std::string desc)
    : ew(ew_), r(r_), omega(omega_), cell(cell_), time_left(tleft), 
      fraction(frac), alive(true), descriptor(desc), random(random_)
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
	descriptor = "killed";
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
	return_value = 1;
    if (desc == "census_born")
	return_value = 2;
    if (desc == "boundary_born")
	return_value = 3;
    if (desc == "vol_emission")
	return_value = 4;
    if (desc == "surface_source")
	return_value = 5;

    // collision event descriptors
    if (desc == "scatter")
	return_value = 100;
    if (desc == "low_weight")
	return_value = 101;
    if (desc == "eff_scatter")
	return_value = 102;
    if (desc == "thom_scatter")
	return_value = 103;
 
    // streaming descriptors
    if (desc == "reflection")
	return_value = 200;
    if (desc == "stream")
	return_value = 201;
    if (desc == "escape")
	return_value = 202;
    if (desc == "cross_boundary")
	return_value = 203;

    // time and census descriptors
    if (desc == "census")
	return_value = 300;

    // death
    if (desc == "killed")
	return_value = 1000;

    // return
    return return_value;
}

//---------------------------------------------------------------------------//
// convert an int into a particle event descriptor

template<class MT>
inline std::string Particle<MT>::get_descriptor(int index)
{
    using std::string;

    // declare return type
    string return_value;

    // born descriptors
    if (index == 1)
	return_value = "born";
    if (index == 2)
	return_value = "census_born";
    if (index == 3)
	return_value = "boundary_born";
    if (index == 4)
	return_value = "vol_emission";
    if (index == 5)
	return_value = "surface_source";

    // collision event descriptors
    if (index == 100)
	return_value = "scatter";
    if (index == 101)
	return_value = "low_weight";
    if (index == 102)
	return_value = "eff_scatter";
    if (index == 103)
	return_value = "thom_scatter";
 
    // streaming descriptors
    if (index == 200)
	return_value = "reflection";
    if (index == 201)
	return_value = "stream";
    if (index == 202)
	return_value = "escape";
    if (index == 203)
	return_value = "cross_boundary";

    // time and census descriptors
    if (index == 300)
	return_value = "census";

    // death
    if (index == 1000)
	return_value = "killed";

    // return
    return return_value;
}

} // end namespace rtt_imc

#endif                          // __imc_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Particle.hh
//---------------------------------------------------------------------------//
