//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Particle_Buffer.hh
 * \author Thomas M. Evans and Mike Buksas
 * \date   Tue May 12 14:34:33 1998
 * \brief  Particle_Buffer and Particle_Stack header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Particle_Buffer_hh
#define rtt_mc_Particle_Buffer_hh

#include "Particle_Stack.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <vector>

#include <mc/config.h>

namespace rtt_mc 
{

// Forward declaration of Particle_IO class.
template<class PT> class Particle_IO;

//===========================================================================//
/*!
 * \class Particle_Buffer
 *
 * \brief Provides a particle buffer that can be used for communication and
 * I/O. 
 *
 * The Particle_Buffer class provides a buffer that stores particles. Its
 * primary use is as a parallel communication buffer.  Particles are added to
 * the buffer using buffer_particle().  Particles are removed from the buffer
 * using add_to_bank.  It is a base class for two derived communication
 * classes: Recv_Particle_Buffer and Send_Particle_Buffer.  The
 * Recv_Particle_Buffer class adds blocking and non-blocking receive
 * functions.  The Send_Particle_Buffer class adds blocking and non-blocking
 * send functions.  All of the particle data is stored in the Particle_Buffer
 * base class.
 *
 * The only services that the Particle Type (PT) need have to use the
 * Particle_Buffer are:
 *
 * \arg std::vector<char> pack();
 * \arg PT(const std::vector<char> &);
 *
 * The first function is a pack function that packs the particle data into a
 * vector<char> array.  The constructor takes a vector<char> array and builds
 * a new particle.
 *
 * The buffer must be told the size of the packed particle (in bytes) with
 * the static set_size_packed_particle() function before buffering functions
 * are used.
 *
 * The default number of particles that can be put into the buffer is 1000.
 * This number can be changed by using the set_maximum_num_particles()
 * function.  The code checks to see that the 32-bit integer limit (number of
 * particles * size of packed particle in bytes) is not exceeded for the
 * total size (in bytes) of the packed particles in the buffer.
 *
 * In general, the derived classes (Recv_Particle_Buffer or
 * Send_Particle_Buffer) will be instantiated directly, even though it is
 * possible to create a Particle_Buffer object.  The base class only provides
 * an interface for adding and removing particles from the buffer; it
 * provides no communication services.  All of the message-passing
 * communication services are in the derived classes.
 *
 * This buffer can also be used by Particle_IO to read/write particles to
 * disk.
 *
 * \sa tstParticle_Buffer.cc for usage examples.
 *
 */
/*!
 * \example mc/test/tstParticle_Buffer.cc
 *
 * Example usage of the Particle_Buffer, Recv_Particle_Buffer, and
 * Send_Particle_Buffer classes.
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
//  7) 21-DEC-01: moved to mc; massive cleanup, removed census particle and
//                comm_buffer
//  8) 10-FEB-03: removed COMPAQ scoping sets because we have found out that
//                this is what the standard dictates; thus we have added base
//                class scoping in the derived classes
//
//===========================================================================//

template<class PT>
class Particle_Buffer
{
    // >>> FRIENDSHIP
    friend class Particle_IO<PT>;

  private:
    // >>> STATIC DATA

    // Maximum number of particles allowed.
    static int max_num_particles;

    // Size of packed particle.
    static int size_packed_particle;

  protected:
    // >>> DATA
    
    // Particle state buffers for receiving.
    std::vector<char> packed_particles;
    
    // Number of particles in the buffer and size of "real" data in the
    // buffer.
    int int_data[2];
    
    // C4_Req communication handles.
    C4::C4_Req comm_int;
    C4::C4_Req comm_char; 

  public:
    // Constructors.
    inline Particle_Buffer();
	
    // Copy Constructor and assignment operators.
    Particle_Buffer(const Particle_Buffer &);
    const Particle_Buffer& operator=(const Particle_Buffer &);

    // Virtual destructor.
    virtual ~Particle_Buffer() { /*...*/ }

    // >>> SET FUNCTIONS

    //! Set the number of particles in the buffer.
    static void set_maximum_num_particles(int);

    //! Set the size of a packed particle.
    static void set_size_packed_particle(int);

    // Empty the Particle_Buffer
    inline void empty_buffer();
    
    // >>> BUFFERING FUNCTIONS

    //! Write a particle to a Particle_Buffer.
    void buffer_particle(const PT &);

    //! Add a Particle_Buffer of particles to a particle bank.
    void add_to_bank(typename Particle_Containers<PT>::Bank &);

    // >>> COMMON PARTICLE COMMUNICATION FUNCTIONS

    // Non-blocking communication services common to sending and receiving.
    void async_free();
    inline bool comm_status() const;

    // >>> ACCESSORS

    //! Get the number of particles stored in the buffer.
    static int get_maximum_num_particles() { return max_num_particles; }

    //! Get the stored size of a packed particle.
    static int get_size_packed_particle() { return size_packed_particle; }

    //! Number of particles in buffer.
    int get_num_particles_in_buffer() const { return int_data[0]; }

    // Check if buffer is empty.
    inline bool is_empty() const; 

    // Check if buffer is full.
    inline bool is_full() const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR PARTICLE BUFFER<PT>::COMM_BUFFER
//---------------------------------------------------------------------------//
/*!
 * \brief Default Particle_Buffer constructor.
 *
 * Initializes an empty Particle_Buffer.
 */
template<class PT>
Particle_Buffer<PT>::Particle_Buffer()  
{
    // set the int_data to zero
    
    // num_particles in buffer
    int_data[0] = 0;

    // size of buffer
    int_data[1] = 0;

    Ensure (packed_particles.size() == int_data[1]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Empty the Particle_Buffer.
 *
 * Sets number of particles in buffer to zero, size of buffer to zero, and
 * empties the buffer.
 *
 * Do not call this function directly after a post_recv.  You must call
 * async_wait (or async_check or async_free) first to make sure the buffer
 * communication has been completed.
 */
template<class PT>
void Particle_Buffer<PT>::empty_buffer()
{
    Insist (!comm_status(), 
	    "Tried to empty an actively communicating buffer.");
 
    // empty the buffer data
    int_data[0] = 0;
    int_data[1] = 0;

    packed_particles.resize(0);

    Ensure (packed_particles.size() == int_data[1]);
    Ensure (int_data[1] == 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check if buffer is empty.
 */
template<class PT>
bool Particle_Buffer<PT>::is_empty() const
{
    // check to see if everything is empty
    return (int_data[0] == 0 && int_data[1] == 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check if buffer is full.
 */
template<class PT>
bool Particle_Buffer<PT>::is_full() const
{
    return (int_data[0] == max_num_particles);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if non-blocking communication is active.
 *
 * The communication buffers are in use until freed by a async_wait,
 * async_check, or async_free.
 */
template<class PT>
bool Particle_Buffer<PT>::comm_status() const
{
    if (!comm_int.inuse())
	return false;
    else if (!comm_char.inuse())
	return false;

    // if we haven't returned then these are still active
    return true;
}

//===========================================================================//
/*!
 * \class Recv_Particle_Buffer
 *
 * This derived class of Particle_Buffer adds blocking and non-blocking
 * receive functions to the Particle_Buffer.
 *
 * The class contains no data.
 */
//===========================================================================//

template<class PT>
class Recv_Particle_Buffer : public Particle_Buffer<PT>
{
  private:
    // Base class typedef for scoping operations.
    typedef Particle_Buffer<PT> Base;

  public:
    // Constructor.
    Recv_Particle_Buffer() {/*...*/}

    // Blocking receive function.
    void recv_buffer(int);

    // Non-blocking receive functions.
    void post_arecv(int);
    void async_wait();
    bool async_check();
};

//===========================================================================//
/*!
 * \class Send_Particle_Buffer
 *
 * This derived class of Particle_Buffer adds blocking and non-blocking
 * send functions to the Particle_Buffer.
 *
 * The class contains no data.
 */
//===========================================================================//

template<class PT>
class Send_Particle_Buffer : public Particle_Buffer<PT>
{
  private:
    // Base class typedef for scoping operations.
    typedef Particle_Buffer<PT> Base;

  public:
    // Constructor.
    Send_Particle_Buffer() {/*...*/}

    // Blocking receive function.
    void send_buffer(int);

    // Non-blocking send functions.
    void post_asend(int);
    void async_wait();
    bool async_check();
};

} // end namespace rtt_mc

#endif                          // rtt_mc_Particle_Buffer_hh

//---------------------------------------------------------------------------//
//                              end of mc/Particle_Buffer.hh
//---------------------------------------------------------------------------//
