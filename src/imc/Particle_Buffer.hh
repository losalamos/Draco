//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Particle_Buffer.hh
 * \author Thomas M. Evans
 * \date   Tue May 12 14:34:33 1998
 * \brief  Particle_Buffer and Particle_Stack header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Particle_Buffer_hh__
#define __imc_Particle_Buffer_hh__

#include "c4/global.hh"
#include "rng/Random.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>

namespace rtt_imc 
{

//===========================================================================//
/*!
 * \class Particle_Stack
 *
 * \brief Stack class for holding particle types (PT).
 *
 * The Particle_Stack class is an implementation of a STL-like stack class.
 * That is, storage is last-in-first-out.  The only difference is that it is
 * defined on std::vector, and it provides and overloaded [] operator for
 * subscripting access.
 *
 * This was originally implemented because of deficiencies in the KCC
 * implementation of std::stack; however, we found that we needed the extra
 * subscripting functionality.
 *
 */
//===========================================================================//

template<class PT>
class Particle_Stack
{
  public:
    // Typedefs.
    typedef typename std::vector<PT>::value_type value_type;
    typedef typename std::vector<PT>::size_type  size_type;

  private:
    // Container holding data in the stack.
    std::vector<PT> c;

  public:
    //! Constructor.
    explicit Particle_Stack(const std::vector<PT> &ct = std::vector<PT>()) 
	: c(ct) {}
    
    //! Query if the stack is empty.
    bool empty() const { return c.empty(); }
    
    //! Return the size of the stack.
    size_type size() const { return c.size(); }
    
    //! Get a reference to the object on the top of the stack.
    value_type& top() { return c.back(); } 

    //! Get a const reference to the object on the top of the stack.
    const value_type& top() const { return c.back(); }
    
    //! Push an object onto the top stack.
    void push(const value_type &x) { c.push_back(x); }

    //! Remove an object from the top of the stack.
    void pop() { c.pop_back(); }

    //! Overloaded operator [] for viewing elements sequentially.
    const value_type& operator[](int i) const { return c[i]; }
};

//===========================================================================//
/*!
 * \class Particle_Buffer
 *
 * \brief Provides particle (PT) communication and persistence services.
 *
 * The Particle_Buffer class provides services for the following:
 *
 * \arg buffering particles for communication
 * \arg sending and receiving particles
 * \arg writing particles to disk
 * \arg bank objects based on the Particle_Stack class
 *
 */
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
//  6) 26-JUL-01: cleaned up file a little, add typedef typename PT::Pack
//                PT_Pack; this typedef is needed because (I'm not 100% sure
//                about the standard) some compilers (g++) have trouble with
//                the syntax PT::Pack::static_function().  This may actually
//                be part of the standard (the gnu developers think so) that
//                KCC is allowing us to get away with.  So remember to use
//                typenames. 
//
//===========================================================================//

template<class PT>
class Particle_Buffer
{
  public:
    // >>> NESTED TYPES

    // Typedefs (this is needed to use static functions of a nested type
    // where the primary type is a template parameter).
    typedef typename PT::Pack    PT_Pack;
    
    // Other typedefs.
    typedef std::vector<double>  sf_double;
    typedef std::vector<int>     sf_int;
    typedef rtt_rng::Sprng       Rnd_Type;
    typedef rtt_rng::Rnd_Control Rnd_Type_Control;
    typedef rtt_dsxx::SP<PT>     SP_PT;
    typedef std::istream         std_istream;
    typedef std::ostream         std_ostream;

    /*!
     * \class Particle_Buffer::Census_Buffer
     * \brief Hold the necessary data for census particles so they can be
     *        written and read from disk.
     */
    struct Census_Buffer
    {
	// Particle state variables.
	sf_double r;
	sf_double omega;
	double    ew;
	double    fraction;
	int       cell;
	Rnd_Type  random;

	// Constructor.
	Census_Buffer(sf_double &, sf_double &, double, double, int,
		      Rnd_Type);

	// Make into a Particle.
	SP_PT make_Particle() const;
    };

    /*!
     * \class Particle_Buffer::Comm_Buffer
     * \brief Native data struct that holds particles for communication.
     */
    struct Comm_Buffer
    {
	// Particle state buffers for receiving.
	double *array_d;
	int    *array_i;
	char   *array_c;

	// C4_Req communication handles.
	C4::C4_Req comm_n;
	C4::C4_Req comm_d;
	C4::C4_Req comm_i;
	C4::C4_Req comm_c;

	// Number of particles in the buffer.
	int n_part;

	// Constructor and destructor.
	inline Comm_Buffer();
	inline ~Comm_Buffer();

	// Copy Constructor and assignment operators.
	inline Comm_Buffer(const Comm_Buffer &);
	inline const Comm_Buffer& operator=(const Comm_Buffer &);
    };

    // Particle banks and containers.
    typedef Particle_Stack<SP_PT> Census;
    typedef Particle_Stack<SP_PT> Bank;

    // Particle communication buffer containers.
    typedef std::vector<Comm_Buffer>    Comm_Vector;
    typedef Particle_Stack<Comm_Buffer> Comm_Bank;

  private:
    // >>> DATA

    // Number of particle doubles that get saved to census (1 less than
    // inflight). 
    int dsize;

    // Number of particle integers.
    int isize;
    
    // Number of particle characters.
    int csize;

    // Number of particles that can be buffered.
    static int buffer_s;

    // Number of particle doubles in the buffer.
    static int buffer_d;

    // Number of particle integers in the buffer.
    static int buffer_i;
    
    // Number of particle char in the buffer.
    static int buffer_c;

  private:
    // >>> Implementation.
    
    // Buffer set functions.
    static void set_buffer(int, int, int);       
    static void set_buffer(int, int, int, int);

  public:
    // Constructors.
    template<class MT>
    Particle_Buffer(const MT &, const Rnd_Type_Control &); 
    Particle_Buffer(int, const Rnd_Type_Control &);
    Particle_Buffer(int, int, int);

    // >>> SET FUNCTIONS

    //! Set the number of particles in the buffer.
    static void set_buffer_size(int);

    // >>> I/O FUNCTIONS.

    //! Write a census particle to disk.
    void write_census(std_ostream &, const PT &) const;

    //! Write a Comm_Buffer of particles to disk.
    void write_census(std_ostream &, Comm_Buffer &) const;

    //! Read a census particle from the disk.
    rtt_dsxx::SP<Census_Buffer> read_census(std_istream &) const;

    // >>> BUFFERING FUNCTIONS

    //! Write a Census_Buffer particle to a Comm_Buffer.
    void buffer_census(Comm_Buffer &, const Census_Buffer &) const;

    //! Write a particle to a Comm_Buffer.
    void buffer_particle(Comm_Buffer &, const PT &) const;

    //! Add a Comm_Buffer of particles to a particle bank.
    void add_to_bank(Comm_Buffer &, Bank &) const;

    // >>> PARTICLE COMMUNICATION FUNCTIONS

    // Blocking send/receives.
    void send_buffer(Comm_Buffer &, int) const;
    rtt_dsxx::SP<Comm_Buffer> recv_buffer(int) const;

    // Non blocking send/receives.
    void asend_buffer(Comm_Buffer &, int) const;
    void post_arecv(Comm_Buffer &, int) const;
    void async_wait(Comm_Buffer &) const;
    bool async_check(Comm_Buffer &) const;
    void async_free(Comm_Buffer &) const;
    bool comm_status(Comm_Buffer &) const;

    // >>> ACCESSORS

    //! Get the number of doubles in a (census) particle.
    int get_dsize() const { return dsize; } 

    //! Get the number of doubles in a census particle.
    int get_census_dsize()   const { return dsize; }

    //! Get the number of doubles in an inflight particle.
    int get_inflight_dsize() const { return dsize + 1; } 

    //! Get the number of integers in a particle.
    int get_isize() const { return isize; }
    
    //! Get the number of characters in a particle.
    int get_csize() const { return csize; }

    //! Get the number of doubles stored in the buffer.
    static int get_buffer_d() { return buffer_d; }

    //! Get the number of integers stored in the buffer.
    static int get_buffer_i() { return buffer_i; }

    //! Get the number of chars stored in the buffer.
    static int get_buffer_c() { return buffer_c; }

    //! Get the number of particles stored in the buffer.
    static int get_buffer_s() { return buffer_s; }
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
