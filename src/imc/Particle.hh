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
// Purpose      : base class which creates and transports particles
//                through a Mesh
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
  // friends and such
    friend class Diagnostic;

public:
  // nested diagnostic class
    class Diagnostic
    {
    private:
      // stream output is sent to
	ostream &output;
    public:
      // constructor
	Diagnostic(ostream &output_) : output(output_) {}

      // diagnostic print functions
	void Print(const Particle<MT> &) const;
	void Print_alive(const Particle<MT> &) const;
	void Print_dead(const Particle<MT> &) const;
    };
	
private:
  // particle weight
    double weight;
  // particle energy
    double energy;
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
    void Stream(double distance)
    {
      // calculate new location
	for (int i = 0; i <= r.size()-1; i++)
	    r[i] = r[i] + distance * omega[i];
    }

  // collision, return a false if particle is absorbed
    bool Collide(const MT &, const Opacity<MT> &);
  // surface crossings, return a false if particle escapes
    bool Surface(const MT &, int);

  // have not yet defined copy constructors or assignment operators
    Particle(const Particle<MT> &);
    Particle<MT>& operator=(const Particle<MT> &);

public:
  // default constructor, explicit to guarantee definition
  // of coord. sys during Particle instantiation
    explicit Particle(const MT &, int, double = 0, double = 0);

  // transport member functions
    void Source(vector<double> &, vector<double> &, const MT &);
    void Transport(const MT &, const Opacity<MT> &, 
		   SP<Diagnostic> = SP<Diagnostic>());

  // other services
    bool Status() const { return alive; }

  // public diagnostic services
    void Print(ostream &) const;
};

// overloaded functions
template<class MT>
ostream& operator<<(ostream &, const Particle<MT> &);

CSPACE

#endif                          // __imctest_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle.hh
//---------------------------------------------------------------------------//
