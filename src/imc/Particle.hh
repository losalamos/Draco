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
#include "ds++/Assert.hh"
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
using std::istream;
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

  // particle buffer class for particle persistence
    class Particle_Buffer
    {
    private:
      // particle data of type double
	double *ddata;
	int dsize;
      // particle data of type int
	int idata;

    public:
      // constructors
	inline Particle_Buffer(const Particle<MT> &);
	inline Particle_Buffer(int, int = 0);

      // destructors
	~Particle_Buffer() { delete [] ddata; }

      // re-calculate particle attributes
	inline vector<double> get_r() const;
	inline vector<double> get_omega() const;
	double get_ew() const { return ddata[0]; }
	double get_frac() const { return ddata[1]; }
	int get_cell() const { return idata; }
	inline Sprng get_rn() const;

      // read and write buffer
	inline void write(ostream &) const;
	inline bool read(istream &);
    };

  // friends and such
    friend class Diagnostic;
    friend class Particle_Buffer;

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
    inline Particle(vector<double>, vector<double>, double, int, Sprng, 
		    double = 1, double = 1);

  // null constructor required as kluge for the STL containers which need a
  // default constructor, this calls an assert(0) so you can't use it
    Particle() { Insist (0, "You tried to default construct a Particle!"); }

  // transport solvers

  // IMC transport step
    void transport(const MT &, const Opacity<MT> &, Tally<MT> &,
		   SP<Diagnostic> = SP<Diagnostic>());

  // other services
    bool status() const { return alive; }
    inline void write_to_census(ostream &) const;

  // public diagnostic services
    void print(ostream &) const;

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// output ascii version of Particle

template<class MT>
inline ostream& operator<<(ostream &output, const Particle<MT> &object)
{
    object.print(output);
    return output;
}

//---------------------------------------------------------------------------//
// dump the Particle_Buffer

template<class MT>
inline ostream& operator<<(ostream &output, const
			   Particle<MT>::Particle_Buffer &object)
{
    object.write(output);
    return output;
}

//---------------------------------------------------------------------------//
// read the Particle_Buffer

template<class MT>
inline bool operator>>(istream &input,  Particle<MT>::Particle_Buffer &object)
{
  // if true we can keep on going, if false we stop
    return object.read(input);
}

//---------------------------------------------------------------------------//
// inline functions for Particle<MT>::Diagnostic
//---------------------------------------------------------------------------//
// Particle<MT>::Diagnostic inline functions

template<class MT>
inline void Particle<MT>::Diagnostic::header() const 
{ 
    output << "*** PARTICLE HISTORY ***" << endl; 
    output << "------------------------" << endl;
}

//---------------------------------------------------------------------------//
// Particle<MT>::Particle_Buffer inline functions
//---------------------------------------------------------------------------//
// constructor for dumping the buffer

template<class MT>
inline Particle<MT>::Particle_Buffer::Particle_Buffer(const Particle<MT>
						      &particle)
{
  // set the size of the dynamic arrays
    dsize = particle.r.size() + 5;
    ddata = new double[dsize];
    
  // assign everything
    int index = 0;
    ddata[index++] = particle.ew;
    ddata[index++] = particle.fraction;
    for (int i = 0; i < particle.omega.size(); i++)
	ddata[index++] = particle.omega[i];
    for (int i = 0; i < particle.r.size(); i++)
	ddata[index++] = particle.r[i];
    Check (index == dsize);
    idata = particle.cell;
}

//---------------------------------------------------------------------------//
// constructor for reading the buffer

template<class MT>
inline Particle<MT>::Particle_Buffer::Particle_Buffer(int dim, int parameter)
{
  // dynamically allocate the size of the double array
    dsize = dim + 5;
    ddata = new double[dsize];
}

//---------------------------------------------------------------------------//
// dump the buffer to an output

template<class MT>
inline void Particle<MT>::Particle_Buffer::write(ostream &output) const
{
  // make sure file exists
    Check (output);

  // dump the output
    output.write(reinterpret_cast<const char *>(ddata), dsize *
		 sizeof(double));
    output.write(reinterpret_cast<const char *>(&idata), sizeof(int));
}

//---------------------------------------------------------------------------//
// get omega from the buffer

template<class MT>
inline vector<double> Particle<MT>::Particle_Buffer::get_omega() const
{
  // make the omega vector
    vector<double> omega;

  // assign the vector
    for (int i = 2; i < 5; i++)
	omega.push_back(ddata[i]);
    
  // return the vector
    return omega;
}

//---------------------------------------------------------------------------//
// get r (position) from the buffer

template<class MT>
inline vector<double> Particle<MT>::Particle_Buffer::get_r() const
{
  // make the r vector
    vector<double> r;

  // assign the vector
    for (int i = 5; i < dsize; i++)
	r.push_back(ddata[i]);

  // return the vector
    return r;
}

//---------------------------------------------------------------------------//
// read the contents

template<class MT>
inline bool Particle<MT>::Particle_Buffer::read(istream &input)
{
  // make sure file exists
    Check (input);
    Check (dsize != 0);

  // read in data
    input.read(reinterpret_cast<char *>(ddata), dsize * sizeof(double));
    if (input.eof())
	return false;
    else
    {
	input.read(reinterpret_cast<char *>(&idata), sizeof(int));
	Check (!input.eof());
    }
    return true;
}

//---------------------------------------------------------------------------//
// Particle<MT> inline functions
//---------------------------------------------------------------------------//
// Particle<MT> constructor

template<class MT>
inline Particle<MT>::Particle(vector<double> r_, vector<double> omega_, 
			      double ew_, int cell_, Sprng random_, 
			      double frac, double tleft)
    : ew(ew_), r(r_), omega(omega_), cell(cell_), time_left(tleft), 
      fraction(frac), alive(true), descriptor("born"), random(random_)
{
  // non-default particle constructor
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

//---------------------------------------------------------------------------//
// write out particle data to an output

template<class MT>
inline void Particle<MT>::write_to_census(ostream &output) const
{
  // make a particle buffer
    Particle_Buffer buffer(*this);

  // dump particle buffer
    buffer.write(output);
}

CSPACE

#endif                          // __imctest_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Particle.hh
//---------------------------------------------------------------------------//
