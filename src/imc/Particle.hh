//----------------------------------*-C++-*----------------------------------//
// Particle.hh
// Thomas M. Evans
// Fri Jan 30 17:04:24 1998
//---------------------------------------------------------------------------//
// @> 
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

#include "Names.hh"
#include "Random.hh"
#include "Opacity.hh"
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
  // random number object
    Random random;
  // stream a distance d
    void stream(double);
  // collision, return a false if particle is absorbed
    bool collide(const MT &, const Opacity<MT> &);
public:
  // default constructor, explicit to guarantee definition
  // of coord. sys during Particle instantiation
    explicit Particle(const MT &, int, double = 0, double = 0);
  // transport member functions
    void source(vector<double> &, const MT &);
    void transport(const MT &, const Opacity<MT> &);
};

CSPACE

#endif                          // __imctest_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle.hh
//---------------------------------------------------------------------------//
