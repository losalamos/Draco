//----------------------------------*-C++-*----------------------------------//
// Particle_Buffer.hh
// Thomas M. Evans
// Tue May 12 14:34:33 1998
//---------------------------------------------------------------------------//
// @> Particle_Buffer class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Particle_Buffer_hh__
#define __imc_Particle_Buffer_hh__

//===========================================================================//
// class Particle_Buffer - 
//
// Purpose : Holds Particles for writing to census or transporting across a
//           cell/processor boundary
//
// revision history:
// -----------------
//  0) original
//  1)  5-26-98 : added temporary Particle_Stack class to account for the
//                deficiency of the KCC 3.3 stack
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/Constants.hh"
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

//===========================================================================//
// class Particle_Stack - 
// Temporary class to account for the KCC 3.3 parser/stack deficiency,
// ie. the KCC 3.3 compiler expects the type to have ==, !=, <= etc defined.
// These constraints should not be placed on the user-defined type.
//===========================================================================//

template<class PT>
class Particle_Stack
{
public:
  // typedefs
    typedef typename list<PT>::value_type value_type;
    typedef typename list<PT>::size_type size_type;

private:
  // container
    list<PT> c;

public:
  // constructor
    explicit Particle_Stack(const list<PT> &ct = list<PT>()) : c(ct) {}
    
  // members
    bool empty() const { return c.empty(); }
    size_type size() const { return c.size(); }
    value_type& top() { return c.back(); } 
    const value_type& top() const { return c.back(); }
    void push(const value_type &x) { c.push_back(x); }
    void pop() { c.pop_back(); }
};

//===========================================================================//
// class Particle_Buffer
//===========================================================================//

template<class PT>
class Particle_Buffer
{
public:
  // abbreviated Particle data from census
    struct Census_Buffer
    {
      // particle state
	vector<double> r;
	vector<double> omega;
	double ew;
	double fraction;
	int cell;
	Sprng random;

      // constructor
	Census_Buffer(vector<double> &, vector<double> &, double, double,
		      int, Sprng);
      // faux default constructor for STL
	Census_Buffer();
    };

  // particle buffer for async receives of particles
    struct Comm_Buffer
    {
      // particle state buffers for receiving
	double array_d[Global::buffer_d];
	int    array_i[Global::buffer_i];
	char   array_c[Global::buffer_c];

      // C4_Req communication handles
	C4_Req comm_n;
	C4_Req comm_d;
	C4_Req comm_i;
	C4_Req comm_c;

      // number of particles in the buffer
	int n_part;
    };

  // standard buffers for particles
    typedef Particle_Stack<PT> Census;
    typedef Particle_Stack<PT> Bank;
    typedef Particle_Stack<PT> Comm_Bank;

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
    void write_census(ostream &, Comm_Buffer &) const;
    SP<Census_Buffer> read_census(istream &) const;

  // fill buffer functions
    void buffer_census(Comm_Buffer &, const Census_Buffer &) const;
    void buffer_particle(Comm_Buffer &, const PT &) const;
  
  // Particle send and receives

  // blocking
    void send_buffer(Comm_Buffer &, int) const;
    SP<Comm_Buffer> recv_buffer(int) const;

  // async
    void asend_bank(Comm_Buffer &, int, Comm_Bank &) const;
    void arecv_bank(Comm_Buffer &, int) const;
    void add_to_bank(Comm_Buffer &, Comm_Bank &) const;

};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR PARTICLE BUFFER
//---------------------------------------------------------------------------//

CSPACE

#endif                          // __imc_Particle_Buffer_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Particle_Buffer.hh
//---------------------------------------------------------------------------//
