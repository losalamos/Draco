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
// Date created : 1-12-97
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
#include <vector>

IMCSPACE

// STL classes used in Particle
using std::vector;

template<class MT>
class Particle
{
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
  // random number object
    Random random;

  // private particle service functions

  // stream a distance d
    void Stream(double);
  // collision, return a false if particle is absorbed
    bool Collide(const MT &, const Opacity<MT> &);

  // private diagnostic functions
    void Print();
    void Print_alive();
    void Print_dead();

public:
  // default constructor, explicit to guarantee definition
  // of coord. sys during Particle instantiation
    explicit Particle(const MT &, int, double = 0, double = 0);

  // transport member functions
    void Source(vector<double> &,vector<double> &, const MT &);
    void Transport(const MT &, const Opacity<MT> &);

  // other services
    bool Status() { return alive; }
};

CSPACE

#endif                          // __imctest_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle.hh
//---------------------------------------------------------------------------//
