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
// 
//===========================================================================//

//===========================================================================//
// struct Particle_Stack - 
//
// Purpose : holds stl types for buffers for Particle types
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Random.hh"
#include "imctest/Opacity.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <stack>
#include <list>

IMCSPACE

// STL classes used in Particle
using std::vector;
using std::string;
using std::ostream;
using std::log;
using std::exp;
using std::stack;
using std::list;

template<class MT> class Particle;

template<class MT>
struct Particle_Stack
{
    typedef stack<Particle<MT>, list<Particle<MT> > > Bank;
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
    template<class T, class C> friend class stack;
    template<class T, class A> friend class list;

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
    Random random;

  // private particle service functions

  // stream a distance d
    inline void stream(double);  

  // stream a distance d
    inline void stream_IMC(const Opacity<MT> &, double);

  // collision, return a false if particle is absorbed
    bool collide(const MT &, const Opacity<MT> &);

  // effective scatter
    void scatter(const MT & );

  // surface crossings, return a false if particle escapes
    bool surface(const MT &, int);

  // null constructor usable only be friends of the class
    Particle();

  // have not yet defined copy constructors or assignment operators
  // Particle(const Particle<MT> &);
  // Particle<MT>& operator=(const Particle<MT> &);

public:
  // explicit constructor
    inline explicit Particle(const MT &, long, double);

  // transport solvers

  // source is temporary until the real source object arrives 
    void source(vector<double> &, vector<double> &, const MT &);

  // IMC transport step
    void transport_IMC(const MT &, const Opacity<MT> &, 
		   SP<Diagnostic> = SP<Diagnostic>());

  // DDMC transport step
    void transport_DDMC(const MT &, const Opacity<MT> &);

  // other services
    bool status() const { return alive; }

  // public diagnostic services
    void print(ostream &) const;
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
    :ew(ew_), r(mesh.get_Coord().get_dim(), 0.0), 
     omega(mesh.get_Coord().get_sdim(), 0.0), cell(0), alive(true), 
     descriptor("born"), random(seed), time_left(0), fraction(1)
{
  // explicit constructor, Particle must be defined with a Mesh
}

template<class MT>
inline void Particle<MT>::stream(double distance)
{
  // calculate new location when Particle streams
    for (int i = 0; i <= r.size()-1; i++)
	r[i] = r[i] + distance * omega[i];
}

template<class MT>
inline void Particle<MT>::stream_IMC(const Opacity<MT> &xs, double distance)
{
  // hardwire minimum energy weight fraction
    double minwt_frac = 0.01;

    double argument = -xs.get_sigeffabs(cell) * distance;
    double min_arg = log(0.1 * minwt_frac);
    if (argument < min_arg) argument = min_arg;

    double factor = exp(argument);
    double new_ew = ew * factor;
    double del_ew = ew - new_ew;

  // Tally::deposit( del_ew, cell );

    fraction *= factor;

    if (fraction < minwt_frac) // kill particle and deposit it energy
    {
      // Tally::deposit( new_ew, cell );
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
