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
// 0) original
// 1) 2-3-98: made Particle a class template parameterized on Mesh_type
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Random.hh"
#include "imctest/Opacity.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>
#include <iostream>

IMCSPACE

// STL classes used in Particle
using std::vector;
using std::string;
using std::ostream;

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
	void print_dist(double, double, int) const;
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
  // status of particle
    bool alive;
  // event type descriptor
    string descriptor;

  // random number object
    Random random;

  // private particle service functions

  // stream a distance d
    inline void stream(double);

  // collision, return a false if particle is absorbed
    bool collide(const MT &, const Opacity<MT> &);

  // surface crossings, return a false if particle escapes
    bool surface(const MT &, int);

  // have not yet defined copy constructors or assignment operators
    Particle(const Particle<MT> &);
    Particle<MT>& operator=(const Particle<MT> &);

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
ostream& operator<<(ostream &, const Particle<MT> &);

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
     descriptor("born"), random(seed)
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

CSPACE

#endif                          // __imctest_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle.hh
//---------------------------------------------------------------------------//
