//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Communicator.hh
 * \author Thomas M. Evans
 * \date   Thu Jul  9 19:16:28 1998
 * \brief  Communicator header file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Communicator_hh__
#define __mc_Communicator_hh__

#include "Particle_Buffer.hh"
#include "Particle_Stack.hh"
#include "ds++/SP.hh"
#include <vector>
#include <algorithm>
#include <iostream>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class Communicator
 *
 * \brief Communicate particles across processor-domain boundaries.
 *
 * Monte Carlo particles are handed to the Communicator class when they enter
 * a cell indexed by a negative integer.  This indicates to the particle type
 * that it has entered a domain-boundary cell.  The absolute value the
 * negative integer is the boundary cell index; that is, a cell that lives on
 * another processor.  At this point the particle can be handed to the
 * communicator.  The communicator stores particles in Particle_Buffer
 * objects for each processor that the local processor communicates with.
 * When these buffers are full the Communicator asynchronously communicates
 * the particles to their waiting processors.  Users have the choice of
 * receiving the buffers and reposting receives or simply receiving the
 * buffers without posting new receives.
 *
 * Communicator instances are built on processor using the
 * rtt_im::Comm_Builder class.  The Communicator class handles both full DD
 * and DD/replication topologies.  It is not required for replication (full
 * replication) topologies.
 *
 * The particle type must support the following services:
 * 
 * \arg int get_cell():  returns the particle's current cell index
 * \arg void set_cell(): set the particle cell index
 */
/*!
 * \example imc/test/tstCommunicator.cc

 * Example usage of the Communicator and Comm_Builder classes.  Note that
 * communicators are only built in non-replication topologies.

 */
// revision history:
// -----------------
//  0) original
//  1) 12-17-98: added print_recv_status() diagnostic to print out the status 
//               of recv_buffer[] 
//  2) 01-03-02: modified to work with the new rtt_mc::Particle_Buffer
//               classes; also moved to the mc package.
// 
//===========================================================================//

template<class PT>
class Communicator
{
  public:
    // Useful typdefs.
    typedef typename Particle_Stack<PT>::Bank      Bank;
    typedef rtt_dsxx::SP<PT>                       SP_PT;
    typedef std::vector<Recv_Particle_Buffer<PT> > sf_Recv_Buffer;
    typedef std::vector<Send_Particle_Buffer<PT> > sf_Send_Buffer;
    typedef std::vector<int>                       sf_int;
    typedef std::vector<std::vector<int> >         vf_int;

  private:
    // Buffers used to receive particles.
    sf_Recv_Buffer recv_buffer;

    // Buffers used to send particle.
    sf_Send_Buffer send_buffer;

    //>>> Data with boundary cell->across processor info.

    // List of nodes this processor receives from and sends to.
    const sf_int recv_nodes;
    const sf_int send_nodes;

    // Temporary last node place holder for boundary cells that live on more
    // than one processor.
    sf_int last_node;

    // List of nodes for each boundary cell.
    const vf_int boundary_node;

    // List of local cells on the respective node for each boundary cell.
    const vf_int boundary_cell;

  private:
    //>>> PRIVATE IMPLEMENTATION.

    // Convert global->local nodes.
    inline int global_to_local(int) const;

  public:
    // Constructor.
    Communicator(const sf_int &, const sf_int &, const vf_int &,
		 const vf_int &);

    // Buffer and send a particle destined for a new processor.
    int communicate(SP_PT);

    // >>> MESSAGE PASSING FUNCTIONS

    // Post receives on all receive buffers.
    void post();

    // Free the buffers.
    void free();

    // Wait on receive buffers followed by a new post.
    int arecv_post(Bank &);

    // Wait on receive buffers.
    int arecv_wait(Bank &);

    // Flush the buffers.
    sf_int flush();
    void flush_all();

    // Send/Recv empty buffers.
    void asend_end();
    void arecv_end();

    // Get the status of send/recv buffers.
    bool asend_status();
    bool arecv_status();

    // >>> ACCESSOR FUNCTIONS

    //! Get the number of nodes to which data can be sent.
    int num_send_nodes() const { return send_nodes.size(); }

    //! Get the number of nodes from which data can be received.
    int num_recv_nodes() const { return recv_nodes.size(); }

    //! Get the number of boundary cells on this processor.
    int num_bound_cells() const { return boundary_node.size(); }

    //! Return a vector of nodes that this processor receives from.
    const sf_int& get_recv_nodes() const { return recv_nodes; }

    //! Return a vector of nodes that this processor sends to.
    const sf_int& get_send_nodes() const { return send_nodes; }

    //! Return the boundary cell - boundary node data.
    const vf_int& get_b_node() const { return boundary_node; }

    //! Return the boundary cell - local cell data.
    const vf_int& get_b_cell() const { return boundary_cell; }

    // Get the number of particles in the send/recv buffers.
    int get_send_size() const;
    int get_recv_size() const;

    // Diagnostic functions.
    void print(std::ostream &) const;
    void arecv_print_status();
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR THE COMMUNICATOR
//---------------------------------------------------------------------------//
/*!
 * \brief Convert a global node index to a local node index.
 */
template<class PT>
int Communicator<PT>::global_to_local(int global_index) const
{
    std::vector<int>::const_iterator itr = std::find(send_nodes.begin(), 
						     send_nodes.end(), 
						     global_index);
    Ensure (itr != send_nodes.end());
    Ensure (*itr == global_index);
    return itr - send_nodes.begin();
}

} // end namespace rtt_mc

#endif                          // __mc_Communicator_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Communicator.hh
//---------------------------------------------------------------------------//
