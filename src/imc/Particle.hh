//----------------------------------*-C++-*----------------------------------//
// Particle.hh
// Thomas M. Evans
// Fri Jan 30 17:04:24 1998
//---------------------------------------------------------------------------//
// @> Particle class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Particle_hh__
#define __imctest_Particle_hh__

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
// 
//===========================================================================//

//===========================================================================//
// struct Particle_Stack - 
//
// Purpose : holds stl types for buffers for Particle types
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Opacity.hh"
#include "imctest/Tally.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <stack>
#include <list>
#include <cassert>

IMCSPACE

// STL classes used in Particle
using std::vector;
using std::string;
using std::ostream;
using std::log;
using std::exp;
using std::stack;
using std::list;

// DRACO classes used in Particle
using RNG::Sprng;

template<class PT>
struct Particle_Stack
{
    typedef stack<PT, list<PT> > Bank;
};

template<class MT>
class Particle
{
public: 
  // nested diagnostic class
    class Diagnostic
    {
    private:
      // stream output is sent to
	ostream &output;
      // boolean for detailed diagnostic
	bool detail;

    public:
      // constructor
	Diagnostic(ostream &output_, bool detail_ = false) 
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

private:
  // particle energy-weight
    double ew;
  // particle location
    vector<double> r;
  // particle direction
    vector<double> omega;
  // particle cell
    int cell;
  // time remaining in this time step
    double time_left;
  // fraction of original energy weight
    double fraction;
  // status of particle
    bool alive;
  // event type descriptor
    string descriptor;

  // random number object
    Sprng random;

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

  // Begin_Doc particle-int.tex
  // Begin_Verbatim 

public:
  // Particle constructor
    inline Particle(const MT &, long, double);

  // null constructor required as kluge for the STL containers which need a
  // default constructor, this calls an assert(0) so you can't use it
    inline Particle();

  // transport solvers

  // IMC transport step
    void transport(const MT &, const Opacity<MT> &, Tally<MT> &,
		   SP<Diagnostic> = SP<Diagnostic>());

  // other services
    bool status() const { return alive; }

  // public diagnostic services
    void print(ostream &) const;

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
inline ostream& operator<<(ostream &output, Particle<MT> &object)
{
    object.print(output);
    return output;
}

//---------------------------------------------------------------------------//
// inline functions for Particle
//---------------------------------------------------------------------------//

// Particle<MT>::Diagnostic inline functions

template<class MT>
inline void Particle<MT>::Diagnostic::header() const 
{ 
    output << "*** PARTICLE HISTORY ***" << endl; 
    output << "------------------------" << endl;
}

// Particle<MT> inline functions

template<class MT>
inline Particle<MT>::Particle(const MT &mesh, long seed, double ew_)
    : ew(ew_), r(mesh.get_Coord().get_dim(), 0.0), 
      omega(mesh.get_Coord().get_sdim(), 0.0), cell(0), alive(true), 
      descriptor("born"), random(seed), time_left(0), fraction(1)
{
  // non-default constructor, Particle must be defined with a Mesh
}

template<class MT>
inline Particle<MT>::Particle()
    : random(-1)
{
  // default constructor for use with stl containers; must provide an
  // initializer for Random because it doesn't have a default constructor 

  // assertion to kill the run if anybody actually tries to use this
    assert (0);
}

template<class MT>
inline void Particle<MT>::stream(double distance)
{
  // calculate new location when Particle streams
    for (int i = 0; i <= r.size()-1; i++)
	r[i] = r[i] + distance * omega[i];
}

template<class MT>
inline void Particle<MT>::stream_IMC(const Opacity<MT> &xs, Tally<MT> &tally,
				     double distance)
{
  // hardwire minimum energy weight fraction
    double minwt_frac = 0.01;

    double argument = -xs.get_sigeffabs(cell) * distance;
    double min_arg = log(0.1 * minwt_frac);
    if (argument < min_arg) argument = min_arg;

    double factor = exp(argument);
    double new_ew = ew * factor;
    double del_ew = ew - new_ew;

    tally.deposit_energy( cell, del_ew );

    fraction *= factor;

    if (fraction < minwt_frac) // kill particle and deposit it energy
    {
	tally.deposit_energy( cell, new_ew );
	descriptor = "killed";
	alive = false;
    }
    else // update particle energy-weight, time_left, and position.
    {
	ew = new_ew;
	time_left -= distance / Global::c;
	for (int i = 0; i <= r.size()-1; i++)
	    r[i] = r[i] + distance * omega[i];
    }
}

CSPACE

#endif                          // __imctest_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle.hh
//---------------------------------------------------------------------------//
