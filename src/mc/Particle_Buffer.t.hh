//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Particle_Buffer.t.hh
 * \author Thomas M. Evans and Mike Buksas
 * \date   Tue May 12 14:34:34 1998
 * \brief  Particle_Buffer member function definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Particle_Buffer.hh"
#include <sstream>

namespace rtt_mc 
{

//===========================================================================//
// PARTICLE_BUFFER BASE CLASS FUNCTIONS
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTORS AND ASSIGNMENT
//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 */
template<class PT>
Particle_Buffer<PT>::Particle_Buffer(const Particle_Buffer &rhs)
    : packed_particles(rhs.packed_particles)
{
    // copy the integer data
    int_data[0] = rhs.int_data[0];
    int_data[1] = rhs.int_data[1];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Overloaded assignment operator.
 */
template<class PT>
const Particle_Buffer<PT>&
Particle_Buffer<PT>::operator=(const Particle_Buffer &rhs)
{
    // check to see if they are the same buffer
    if (&rhs == this)
	return *this;

    // assign the vector
    packed_particles = rhs.packed_particles;

    // copy the integer data
    int_data[0] = rhs.int_data[0];
    int_data[1] = rhs.int_data[1];

    return *this;
}

//---------------------------------------------------------------------------//
// STATIC PRIVATE DATA MEMBERS
//---------------------------------------------------------------------------//

//! Maximum number of particles allowed in a Particle_Buffer.
template<class PT> int Particle_Buffer<PT>::max_num_particles    = 1000;

//! Size of packed particles.
template<class PT> int Particle_Buffer<PT>::size_packed_particle = 0;

//---------------------------------------------------------------------------//
// STATIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

template<class PT>
void Particle_Buffer<PT>::set_maximum_num_particles(int num)
{
    Require (num > 0);
    Require (size_packed_particle > 0);

    // eight byte int type is determined by configure
    EIGHT_BYTE_INT_TYPE max_bytes = static_cast<EIGHT_BYTE_INT_TYPE>(num) * 
	size_packed_particle;

    EIGHT_BYTE_INT_TYPE limit     = static_cast<EIGHT_BYTE_INT_TYPE>(
	2147483648U);

    Insist (max_bytes < limit, "Buffer size exceeds 32-bit int limit.");

    max_num_particles = num;
}

//---------------------------------------------------------------------------//

template<class PT>
void Particle_Buffer<PT>::set_size_packed_particle(int size)
{
    Require (size > 0);

    size_packed_particle = size;
}

//---------------------------------------------------------------------------//
// BUFFERING FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Add a particle to the Particle_Buffer.
 *
 * The particle type (PT) must be capable of packing itself into a
 * vector<char> (through the function PT::pack()).
 *
 * The Particle_Buffer communicators must have finished any operations before
 * this function can be invoked.  An assertion is thrown if the communicator
 * requests are active.
 */
template<class PT>
void Particle_Buffer<PT>::buffer_particle(const PT &particle)
{
    using std::vector;

    Insist (!comm_status(), 
	    "Tried to buffer a particle in an actively communicating buffer.");

    // check to make sure the buffer isn't full
    Insist (int_data[0] < max_num_particles, 
	    "Tried to overfill Particle_Buffer.");

    Check (int_data[0] >= 0);
    Check (int_data[0] * size_packed_particle == int_data[1]);

    // pack up the particle
    vector<char> pack_PT = particle.pack();
    Check (pack_PT.size() == size_packed_particle);
    
    // if the size of the buffer is equal to int_data[1], then we can simply
    // insert the packed particle onto the end of the buffer
    if (packed_particles.size() == int_data[1])
    {
	// add particle to the end of the buffer
	packed_particles.insert(packed_particles.end(), 
				pack_PT.begin(), pack_PT.end());
	
	// update the number of particles and size of buffer
	int_data[1] += size_packed_particle;
	Check (packed_particles.size() == int_data[1]);

	Check (packed_particles.size() <=    
	       max_num_particles * size_packed_particle);
    }

    // if the size of the buffer is equal to the maximum then we have to copy
    // data into the buffer, this would be the case after a post_arecv
    else if (packed_particles.size() == 
	     max_num_particles * size_packed_particle)
    {
	// copy data into the buffer, the actual amount of data in the
	// buffer is determined by int_data[1]
	vector<char>::iterator itr = packed_particles.begin() + int_data[1];
	vector<char>::iterator end = itr + size_packed_particle;
	for (vector<char>::iterator p = pack_PT.begin(); itr != end; itr++)
	    *itr = *p++; 
	
	// move the buffer size indicator forward
	int_data[1] += size_packed_particle;

	Check (int_data[1] <= packed_particles.size());
    }

    // if the above conditions are untrue than something is probably hosed,
    // most likely the buffer has been tampered with in the middle of
    // send/recv calling
    else
    {
	std::ostringstream message;
	message << "Inconsistency in Particle_Buffer. "
		<< "Possibilities are: calling set_maximum_num_particles "
		<< "in the middle of async send/recv or calling empty_buffer "
		<< "in the middle of asnyc send/recv.";
	throw rtt_dsxx::assertion(message.str().c_str());
    }

    // update the number of particles in the buffer
    int_data[0]++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack particles from buffer and add them to a bank.
 *
 * The Particle_Buffer communicators must have finished any operations before
 * this function can be invoked.  An assertion is thrown if the communicator
 * requests are active.
 *
 * The buffer is emptied after the particles are added to the bank.
 *
 * The particle type (PT) must be capable of unpacking itself from a
 * vector<char> with a constructor (PT::PT(vector<char>)).
 *
 * \param bank Particle_Stack<PT>::Bank, this is simple a typedef to
 * Particle_Stack<SP<PT>> 
 */
template<class PT>
void Particle_Buffer<PT>::add_to_bank(
    typename Particle_Containers<PT>::Bank &bank)
{
    using rtt_dsxx::SP;
    using std::vector;

    Insist (!comm_status(), 
	    "Tried to empty an actively communicating buffer.");

    Insist (size_packed_particle > 0, "Undefined packed particle size.");
    
    Require (int_data[0] >= 0 && int_data[0] <= max_num_particles);
    Require (int_data[0] * size_packed_particle == int_data[1]);
    Require (int_data[1] <= packed_particles.size());

    // iterator to packed_particles
    vector<char>::const_iterator itr = packed_particles.begin();

    // a packed particle
    vector<char> packed_PT(size_packed_particle);
    SP<PT>       particle;

    // number of particles in bank initially
    int bank_size = bank.size();

    // loop through and add the Particles to the Bank.
    for (int n = 0; n < int_data[0]; n++)
    {
	// fill up a packed particle
	for (int i = 0; i < size_packed_particle; i++, itr++)
	    packed_PT[i] = *itr;

	// build the new particle
	particle = new PT(packed_PT);
	Check (particle);

	// add it to the bank
	bank.push(particle);
    }
    Ensure (bank_size + int_data[0] == bank.size());
    Ensure (packed_particles.size() == int_data[1] ? 
	    itr == packed_particles.end() : true);
    
    // empty the buffer
    empty_buffer();

    Ensure (is_empty());
}

//---------------------------------------------------------------------------//
// COMMON COMMUNICATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Free the communication buffers.
 *
 * This function cancels any outstanding communications and should be used
 * with extreme care.  The state of the buffer is effectively undefined after
 * calling async_free.
 */
template<class PT>
void Particle_Buffer<PT>::async_free()
{
    // free the C4_Req objects from Async posts, note that these must be
    // reassigned with a post request before they can be received or tested 
    comm_int.free();
    comm_char.free();

    Ensure (!comm_status());
}

//===========================================================================//
// RECV_PARTICLE_BUFFER DERIVED CLASS FUNCTIONS
//===========================================================================//
//---------------------------------------------------------------------------//
// BLOCKING PARTICLE COMMUNICATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Receive particles into the buffer.
 *
 * This function receives particles from proc using \b blocking communication
 * into the Recv_Particle_Buffer.  The Recv_Particle_Buffer \b must be empty
 * (this condition is checked and an rtt_dsxx::assertion is thrown if false).
 */
template<class PT>
void Recv_Particle_Buffer<PT>::recv_buffer(int proc)
{
    using C4::Recv;

    // buffer must be empty
    Insist (is_empty(), "Tried to receive particles in a non-empty buffer.");

    // receive the ints
    Recv <int> (&int_data[0], 2, proc, 100);
    Check (int_data[0] >= 0);
    Check (int_data[1] >= 0);

    // resize the buffer
    packed_particles.resize(int_data[1]);

    // receive the buffer
    Recv <char> (&packed_particles[0], int_data[1], proc, 101);
}
    
//---------------------------------------------------------------------------//
// NON-BLOCKING PARTICLE COMMUNICATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Post non-blocking receives.
 *
 * Posts non-blocking receives for the particles in the Particle_Buffer.  The
 * Particle_Buffer must be empty; if not an rtt_dsxx::assertion is thrown.
 *
 * Before using the Particle_Buffer for any operations after a post_arecv it
 * is required to check that the buffer has been received by using
 * async_check or async_wait.
 */
template<class PT>
void Recv_Particle_Buffer<PT>::post_arecv(int proc)
{
    using C4::RecvAsync;

    // buffer must be empty
    Insist (is_empty(), "Tried to post receive in a non-empty buffer.");

    // resize the buffer to max size, for async receives the buffer size is
    // set to max, the actual data copied into the buffer on a receive is
    // given by int_data[1]
    packed_particles.resize(get_size_packed_particle() * 
			    get_maximum_num_particles());

    // post receives for int data
    RecvAsync(comm_int, &int_data[0], 2, proc, 100);
    
    // post receives for the buffer
    RecvAsync(comm_char, &packed_particles[0], packed_particles.size(), 
	      proc, 101);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait on a non-blocking receive operation.
 * 
 * This function waits on a non-blocking receive.  This function, or
 * async_check, must be called before doing work with the Particle_Buffer
 * after a post_arecv.
 */
template<class PT> 
void Recv_Particle_Buffer<PT>::async_wait()
{
    // wait on receive buffers to make sure they are full
    comm_int.wait();
    comm_char.wait();

    Ensure (int_data[1] <= packed_particles.size());
    Ensure (int_data[0] >= 0 && int_data[0] <= get_maximum_num_particles());
    Ensure (int_data[0] * get_size_packed_particle() == int_data[1]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check to see if a non-blocking receive is complete.
 *
 * This function checks on a non-blocking receive.  This function, or
 * async_wait, should be called before doing work with the Particle_Buffer
 * after a post_arecv.
 */
template<class PT>
bool Recv_Particle_Buffer<PT>::async_check() 
{
    using std::vector;

    // tag to check what is in; we want all or nothing
    vector<int> arrived(4, 0);
    int total = 0;
    int count = 0;
    
    // check to see if the buffers have been received
    do
    {
	// check comm_int
	if (arrived[0] == 0)
	    if (comm_int.complete())
	    {
		arrived[0] = 1;
		total++;
	    }
	
	// check comm_char
	if (arrived[1] == 0)
	    if (comm_char.complete())
	    {
		arrived[1] = 1;
		total++;
	    }

	// accumulate count
	count++;
	
    } while (total > 0 && total < 2);

    Ensure (total == 0 || total == 2);
    return total;
}

//===========================================================================//
// SEND_PARTICLE_BUFFER DERIVED CLASS FUNCTIONS
//===========================================================================//
//---------------------------------------------------------------------------//
// BLOCKING PARTICLE COMMUNICATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a blocking send of the particles in the buffer.
 *
 * This sends the particles stored in the Particle_Buffer to processor
 * indicated by proc using \b blocking communication.  The buffer is emptied
 * after the send.
 */
template<class PT>
void Send_Particle_Buffer<PT>::send_buffer(int proc)
{
    using C4::Send;

    // send out the integer data
    Send <int>(&int_data[0], 2, proc, 100);
    
    // send out the buffer
    Send <char>(&packed_particles[0], int_data[1], proc, 101);

    // empty the Particle_Buffer
    empty_buffer();

    Ensure (is_empty());
}

//---------------------------------------------------------------------------//
// NON-BLOCKING PARTICLE COMMUNICATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Post non-blocking send for the Particle_Buffer.
 *
 * This sends the packed particles to proc using non-blocking communication.
 * It does not empty the buffer after sending.  However, the buffer is
 * automatically emptied after calling async_wait() or async_check()
 * (provided the check returns true).
 */
template<class PT>
void Send_Particle_Buffer<PT>::post_asend(int proc)
{
    using C4::SendAsync;

    // send out the integer data
    SendAsync(comm_int, &int_data[0], 2, proc, 100);
    
    // send out the buffer
    SendAsync(comm_char, &packed_particles[0], int_data[1], proc, 101);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait on a non-blocking send operation.
 * 
 * This function waits on a non-blocking send.  This function, or
 * async_check, must be called before doing work with the Particle_Buffer
 * after a post_asend.  The buffer is emptied after the wait is completed.
 */
template<class PT> 
void Send_Particle_Buffer<PT>::async_wait()
{
    // wait on receive buffers to make sure they are full
    comm_int.wait();
    comm_char.wait();

    // empty the buffer
    empty_buffer();

    Ensure (is_empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check to see if a non-blocking send is complete.
 *
 * This function checks on a non-blocking send.  This function, or
 * async_wait, should be called before doing work with the Particle_Buffer
 * after a post_asend.  If the send communication is complete, the buffer is
 * emptied.  If the send is not complete the buffer is untouched.
 */
template<class PT>
bool Send_Particle_Buffer<PT>::async_check() 
{
    using std::vector;

    // tag to check what is in; we want all or nothing
    vector<int> arrived(4, 0);
    int total = 0;
    int count = 0;
    
    // check to see if the buffers have been received
    do
    {
	// check comm_int
	if (arrived[0] == 0)
	    if (comm_int.complete())
	    {
		arrived[0] = 1;
		total++;
	    }
	
	// check comm_char
	if (arrived[1] == 0)
	    if (comm_char.complete())
	    {
		arrived[1] = 1;
		total++;
	    }

	// accumulate count
	count++;
	
    } while (total > 0 && total < 2);

    Ensure (total == 0 || total == 2);

    // if the operation is complete, empty the buffer
    if (total)
	empty_buffer();

    return total;
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Particle_Buffer.t.hh
//---------------------------------------------------------------------------//
