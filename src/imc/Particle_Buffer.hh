//----------------------------------*-C++-*----------------------------------//
// Particle_Buffer.hh
// Thomas M. Evans
// Tue May 12 14:34:33 1998
//---------------------------------------------------------------------------//
// @> Particle_Buffer class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Particle_Buffer_hh__
#define __imctest_Particle_Buffer_hh__

//===========================================================================//
// class Particle_Buffer - 
//
// Purpose : Holds Particles for writing to census or transporting across a
//           cell/processor boundary
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "rng/Random.hh"
#include "ds++/Assert.hh"

IMCSPACE

using RNG::Sprng;

template<class MT>
class Particle_Buffer
{
public:
  // standard buffers for particles
    typedef stack<PT, list<PT> > census;
    typedef stack<PT, list<PT> > bank;
    typedef stack<PT, list<PT> > comm_bank;

  // abbreviated Particle data from census
    struct Census_Particle
    {
      // particle state
	vector<double> r;
	vector<double> omega;
	double ew;
	double fraction;
	int cell;
	Sprng random;

      // constructor
	Census_Particle(vector<double> &, vector<double> &, double, double,
			int, Sprng)
    };

private:
  // data of type double size (number of elements)
    int dsize;
  // data of type int size (number of elements)
    int isize;

public:
  // constructor
    template<class MT> explicit Particle_Buffer(const MT &);

  // io functions
    inline void write_census(const PT &);
    SP<Census_Particle> read_census();
};

CSPACE

#endif                          // __imctest_Particle_Buffer_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle_Buffer.hh
//---------------------------------------------------------------------------//
