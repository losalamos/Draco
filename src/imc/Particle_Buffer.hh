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
#include "imctest/Constants.hh"
#include "c4/global.hh"
#include "rng/Random.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <stack>
#include <list>

IMCSPACE

// draco components
using RNG::Sprng;
using RNG::Rnd_Control;
using C4::node;
using C4::nodes;
using C4::C4_Req;
using C4::SendAsync;
using C4::RecvAsync;

// STL components
using std::vector;
using std::list;
using std::stack;
using std::istream;
using std::ostream;

template<class PT>
class Particle_Buffer
{
public:
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
			int, Sprng);
      // faux default constructor for STL
	Census_Particle();
    };

  // particle buffer for async receives of particles
    struct Comm_Buffer
    {
      // particle state buffers
	double array_d[Global::buffer_d];
	int array_i[Global::buffer_i];
	char array_c[Global::buffer_c];

      // number of particles in the buffer
	int n_part;
    };

  // standard buffers for particles
    typedef stack<PT, list<PT> > Census;
    typedef stack<PT, list<PT> > Bank;
    typedef stack<PT, list<PT> > Comm_bank;

private:
  // data of type double size (number of elements)
    int dsize;
  // data of type int size (number of elements)
    int isize;
  // data of type char size (number of bytes of random number state)
    int csize;

public:
  // constructor
    template<class MT> Particle_Buffer(const MT &, const Rnd_Control &);

  // io functions
    void write_census(ostream &, const PT &) const;
    void write_census(const Census &) const;
    SP<Census_Particle> read_census(istream &);

  // Particle send and receives
    void send_bank(C4_Req &, int, Comm_bank &) const;
    SP<Comm_Buffer> recv_bank(C4_Req &, int) const;
    void add_to_bank(SP<Comm_Buffer>, Comm_bank &) const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR PARTICLE BUFFER
//---------------------------------------------------------------------------//

CSPACE

#endif                          // __imctest_Particle_Buffer_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle_Buffer.hh
//---------------------------------------------------------------------------//
