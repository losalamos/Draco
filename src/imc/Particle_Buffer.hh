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
//  2)  6-11-98 : moved the buffer sizes into Particle_Buffer as private:
//                static variables, they can only be accessed by accessor
//                functions from the outside
//  3)  7-30-98 : add free_arecv function to free Comm_Buffers that are
//                waiting on an async receive
//  4)  4-30-99 : fixed set_buffer_size() function, it had an integer
//                division error that hit us when we reduced the buffer size 
//  5) 7-MAR-00 : added Census_Buffer::make_Particle() function that converts 
//                a Census_Buffer object into a particle.
//
//===========================================================================//

#include "c4/global.hh"
#include "rng/Random.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>

namespace rtt_imc 
{

//===========================================================================//
// class Particle_Stack - 
// Temporary class to account for the KCC 3.3 parser/stack deficiency,
// ie. the KCC 3.3 compiler expects the type to have ==, !=, <= etc defined.
// These constraints should not be placed on the user-defined type. 
// NOTE: The Particle_Stack class is now our preferred class because we have
// added subscripting for sequential access for looking at the data
//===========================================================================//

template<class PT>
class Particle_Stack
{
  public:
    // typedefs
    typedef typename std::vector<PT>::value_type value_type;
    typedef typename std::vector<PT>::size_type size_type;

  private:
    // container
    std::vector<PT> c;

  public:
    // constructor
    explicit Particle_Stack(const std::vector<PT> &ct = std::vector<PT>()) 
	: c(ct) {}
    
    // members
    bool empty() const { return c.empty(); }
    size_type size() const { return c.size(); }
    value_type& top() { return c.back(); } 
    const value_type& top() const { return c.back(); }
    void push(const value_type &x) { c.push_back(x); }
    void pop() { c.pop_back(); }

    // overloaded operator () for viewing elements sequentially
    const value_type& operator[](int i) const { return c[i]; }
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
	std::vector<double> r;
	std::vector<double> omega;
	double ew;
	double fraction;
	int cell;
	rtt_rng::Sprng random;

	// constructor
	Census_Buffer(std::vector<double> &, std::vector<double> &, double,
		      double, int, rtt_rng::Sprng);
	// faux default constructor for STL
	Census_Buffer();

	// make into a Particle
	rtt_dsxx::SP<PT> make_Particle() const;
    };

    // particle buffer for async receives of particles
    struct Comm_Buffer
    {
	// particle state buffers for receiving
	double *array_d;
	int    *array_i;
	char   *array_c;

	// C4_Req communication handles
	C4::C4_Req comm_n;
	C4::C4_Req comm_d;
	C4::C4_Req comm_i;
	C4::C4_Req comm_c;

	// number of particles in the buffer
	int n_part;

	// inline default constructor
	inline Comm_Buffer();
	inline ~Comm_Buffer();

	// inline Copy Constructor and assignment operators
	inline Comm_Buffer(const Comm_Buffer &);
	inline const Comm_Buffer& operator=(const Comm_Buffer &);
    };

    // standard buffers for particles
    typedef Particle_Stack<rtt_dsxx::SP<PT> > Census;
    typedef Particle_Stack<rtt_dsxx::SP<PT> > Bank;

    // standard buffers for Comm_Buffers
    typedef std::vector<Comm_Buffer> Comm_Vector;
    typedef Particle_Stack<Comm_Buffer> Comm_Bank;

  private:
    // data of type double size (number of elements) saved to census
    int dsize;
    // data of type int size (number of elements)
    int isize;
    // data of type char size (number of bytes of random number state)
    int csize;

    // static buffer sizes
    static int buffer_s;
    static int buffer_d;
    static int buffer_i;
    static int buffer_c;

  public:
    // constructor
    template<class MT>
    Particle_Buffer(const MT &, const rtt_rng::Rnd_Control &);  
    Particle_Buffer(int, int, int);

    // buffer sizing and accessor functions
    static void set_buffer(int, int, int);
    static void set_buffer(int, int, int, int);
    static void set_buffer_size(int);
    static int get_buffer_d() { return buffer_d; }
    static int get_buffer_i() { return buffer_i; }
    static int get_buffer_c() { return buffer_c; }
    static int get_buffer_s() { return buffer_s; }

    // io functions
    void write_census(std::ostream &, const PT &) const;
    void write_census(std::ostream &, Comm_Buffer &) const;
    rtt_dsxx::SP<Census_Buffer> read_census(std::istream &) const;

    // fill and get buffer functions
    void buffer_census(Comm_Buffer &, const Census_Buffer &) const;
    void buffer_particle(Comm_Buffer &, const PT &) const;
    void add_to_bank(Comm_Buffer &, Bank &) const;

    // Particle send and receives

    // blocking
    void send_buffer(Comm_Buffer &, int) const;
    rtt_dsxx::SP<Comm_Buffer> recv_buffer(int) const;

    // async
    void asend_buffer(Comm_Buffer &, int) const;
    void post_arecv(Comm_Buffer &, int) const;
    void async_wait(Comm_Buffer &) const;
    bool async_check(Comm_Buffer &) const;
    void async_free(Comm_Buffer &) const;
    bool comm_status(Comm_Buffer &) const;

    // accessor functions
    int get_dsize() const { return dsize; } 
    int get_isize() const { return isize; }
    int get_csize() const { return csize; }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR PARTICLE BUFFER<PT>::Comm_Buffer
//---------------------------------------------------------------------------//
// default constructor

template<class PT>
inline Particle_Buffer<PT>::Comm_Buffer::Comm_Buffer()
    : array_d(new double[get_buffer_d()]), 
      array_i(new int[get_buffer_i()]),
      array_c(new char[get_buffer_c()]),
      n_part(0)
{
    // dynamically sized values to the Globally determined buffer sizes
}

//---------------------------------------------------------------------------//
// destructor to reclaim our memory

template<class PT>
inline Particle_Buffer<PT>::Comm_Buffer::~Comm_Buffer()
{
    // get back our dynamically allocated memory
    delete [] array_d;
    delete [] array_i;
    delete [] array_c;
}

//---------------------------------------------------------------------------//
// copy constructor

template<class PT>
inline Particle_Buffer<PT>::Comm_Buffer::Comm_Buffer(const Comm_Buffer &rhs)
    : array_d(new double[get_buffer_d()]), 
      array_i(new int[get_buffer_i()]),
      array_c(new char[get_buffer_c()]),
      n_part(rhs.n_part)
{
    // copy double data
    for (int i = 0; i < get_buffer_d(); i++)
	array_d[i] = rhs.array_d[i];

    // copy int data
    for (int i = 0; i < get_buffer_i(); i++)
	array_i[i] = rhs.array_i[i];

    // copy char data
    for (int i = 0; i < get_buffer_c(); i++)
	array_c[i] = rhs.array_c[i];
}

//---------------------------------------------------------------------------//
// overloaded assignment operator

template<class PT>
inline const typename Particle_Buffer<PT>::Comm_Buffer&
Particle_Buffer<PT>::Comm_Buffer::operator=(const Comm_Buffer &rhs)
{
    // check to see if they are the same buffer
    if (&rhs == this)
	return *this;

    // we know that these objects are the same size, just do the assignment
    n_part = rhs.n_part;

    // copy double data
    for (int i = 0; i < get_buffer_d(); i++)
	array_d[i] = rhs.array_d[i];

    // copy int data
    for (int i = 0; i < get_buffer_i(); i++)
	array_i[i] = rhs.array_i[i];

    // copy char data
    for (int i = 0; i < get_buffer_c(); i++)
	array_c[i] = rhs.array_c[i];

    // return for concatenated calls
    return *this;
}

} // end namespace rtt_imc

#endif                          // __imc_Particle_Buffer_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Particle_Buffer.hh
//---------------------------------------------------------------------------//
