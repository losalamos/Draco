//----------------------------------*-C++-*----------------------------------//
// Particle.hh 
// Thomas M. Evans
// Fri Jan 30 17:04:24 1998
//---------------------------------------------------------------------------//
// @> Particle class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Particle_hh__
#define __imc_Particle_hh__

//===========================================================================//
// class Particle - 
//
// Purpose : base class which creates and transports particles
//           through a Mesh
//
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
// 
//===========================================================================//

//===========================================================================//
// struct Particle_Stack - 
//
// Purpose : holds stl types for buffers for Particle types
// 
//===========================================================================//

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

template<class MT>
class Particle
{
  public: 
    // nested diagnostic class
    class Diagnostic
    {
      private:
	// stream output is sent to
	std::ostream &output;
	// boolean for detailed diagnostic
	bool detail;

      public:
	// constructor
	Diagnostic(std::ostream &output_, bool detail_ = false) 
	    : output(output_), detail(detail_) {}

	// switches
	bool detail_status() const { return detail; }

	// diagnostic print functions
	void print(const Particle<MT> &) const;
	void print_alive(const Particle<MT> &) const;
	void print_dead(const Particle<MT> &) const;
	void print_dist(double, double, double, int) const;
	void print_xs(const Opacity<MT> &, int) const;

	// inline output formatters
	inline void header() const;
    };

    // friends and such
    friend class Diagnostic;
    template<class PT> friend class Particle_Buffer;

  private:
    // particle energy-weight
    double ew;
    // particle location
    std::vector<double> r;
    // particle direction
    std::vector<double> omega;
    // particle cell
    int cell;
    // time remaining in this time step
    double time_left;
    // fraction of original energy weight
    double fraction;
    // status of particle
    bool alive;
    // event type descriptor
    std::string descriptor;

    // random number object
    rtt_rng::Sprng random;

    // private particle service functions

    // stream a distance d
    inline void stream(double);  

    // stream a distance d
    inline void stream_IMC(const Opacity<MT> &, Tally<MT> &, double);

    // collision, return a false if particle is absorbed
    bool collide(const MT &, const Opacity<MT> &);

    // effective scatter
    void scatter(const MT & );

    // surface crossings, return a false if particle escapes
    bool surface(const MT &, int);

    // IMC transport step
    void trans_IMC(const MT &, const Opacity<MT> &);
    
    // DDMC transport step
    void trans_DDMC(const MT &, const Opacity<MT> &);

  public:
    // Particle constructor
    inline Particle(std::vector<double>, std::vector<double>, double, int,
		    rtt_rng::Sprng, double = 1, double = 1, 
		    std::string = "born");

    // null constructor required as kluge for the STL containers which need a
    // default constructor, this calls an assert(0) so you can't use it
    Particle() { Insist (0, "You tried to default construct a Particle!"); }

    // transport solvers

    // IMC transport step
    void transport(const MT &, const Opacity<MT> &, Tally<MT> &,
		   rtt_dsxx::SP<Diagnostic> = rtt_dsxx::SP<Diagnostic>()); 

    // other services

    // particle status functions
    bool status() const { return alive; }
    void reset_status() { alive = true; }
    void kill_particle() { alive = false; }

    // accessors
    int get_cell() const { return cell; }
    double get_ew() const { return ew; }
    const rtt_rng::Sprng& get_random() const { return random; }

    // set functions for sourcing particles
    void set_random(rtt_rng::Sprng &ran) { random = ran; }
    void set_time_left(double t) { time_left = t; }
    void set_descriptor(std::string s) { descriptor = s; }
    void set_ew(double new_ew) { ew = new_ew; }
    void set_cell(int new_cell) { cell = new_cell; }

    // transport descriptors
    std::string desc() const { return descriptor; }
    inline static int get_index(std::string);
    inline static std::string get_descriptor(int);

    // public diagnostic services
    void print(std::ostream &) const;

    // overloaded operators
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
    // hardwire minimum energy weight fraction
    double minwt_frac = 0.01;
    double min_arg    = log(0.1 * minwt_frac);

    double argument = -xs.get_sigeffabs(cell) * distance;
    if (argument < min_arg) argument = min_arg;

    double factor = exp(argument);
    double new_ew = ew * factor;
    double del_ew = ew - new_ew;

    tally.deposit_energy( cell, del_ew );
    if ( xs.get_sigeffabs(cell) != 0)
	tally.accumulate_ewpl( cell, del_ew / xs.get_sigeffabs(cell) );

    fraction *= factor;

    if (fraction < minwt_frac) // kill particle and deposit it energy
    {
	tally.deposit_energy( cell, new_ew );
	tally.accum_n_killed();
	tally.accum_ew_killed( new_ew );
	descriptor = "killed";
	alive = false;
    }
    else // update particle energy-weight, time_left, and position.
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
