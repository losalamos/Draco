//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Communicator.hh
 * \author Thomas M. Evans
 * \date   Thu Jul  9 19:16:28 1998
 * \brief  Communicator header file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Communicator_hh__
#define __imc_Communicator_hh__

#include "Particle_Buffer.hh"
#include "ds++/SP.hh"
#include <vector>
#include <algorithm>
#include <iostream>

namespace rtt_imc 
{

//===========================================================================//
/*!
 * \class Communicator
 
 * \brief Communicate particles across processor-domain boundaries.

 * IMC particles are handed to the Communicator class when they enter a cell
 * indexed by a negative integer.  This confirms to the particle (represented
 * by the rtt_imc::Particle class) that it has entered a domain-boundary
 * cell.  At this point the particle can be handed to the communicator.  The
 * communicator stores particles in particle buffers (using the
 * rtt_imc::Particle_Buffer class) for each processor that the local
 * processor communicates with.  When these buffers are full the Communicator
 * asynchronously communicates the particles to their waiting processors.
 * Users have the choice of receiving the buffers are reposting receives or
 * simply receiving the buffers without posting new receives.

 * Communicator instances are built on processor using the
 * rtt_imc::Comm_Builder class.  The Communicator class handles both full DD
 * and DD/replication topologies.  It is not required for replication (full
 * replication) topologies.

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
// 
//===========================================================================//

template<class PT>
class Communicator
{
  public:
    // Useful typdefs.
    typedef typename Particle_Buffer<PT>::Bank Bank;
    typedef rtt_dsxx::SP<PT>                   SP_PT;
    typedef std::vector<int>                   sf_int;
    typedef std::vector<std::vector<int> >     vf_int;

  private:
    // Comm_Buffers for particle tranport.
    typename Particle_Buffer<PT>::Comm_Vector recv_buffer;
    typename Particle_Buffer<PT>::Comm_Vector send_buffer;

    //>>> Data with boundary cell->across processor info.

    // List of nodes this processor receives from and sends to.
    sf_int recv_nodes;
    sf_int send_nodes;

    // Temporary last node place holder for boundary cells that live on more
    // than one processor.
    sf_int last_node;

    // List of nodes for each boundary cell.
    vf_int boundary_node;

    // List of local cells on the respective node for each boundary cell.
    vf_int boundary_cell;

  private:
    //>>> PRIVATE IMPLEMENTATION.

    // Convert global->local nodes.
    inline int global_to_local(int) const;

  public:
    // Constructor.
    explicit Communicator(const sf_int &, const sf_int &, const vf_int &,
			  const vf_int &);

    // Buffer a particle destined for a new processor.
    int communicate(const Particle_Buffer<PT> &, SP_PT);

    // Bessage passing (async and sync) operations.
    void post(const Particle_Buffer<PT> &);
    sf_int flush(const Particle_Buffer<PT> &);
    void flush_all(const Particle_Buffer<PT> &);
    void free(const Particle_Buffer<PT> &);
    int arecv_post(const Particle_Buffer<PT> &, Bank &);
    void arecv_wait(const Particle_Buffer<PT> &, Bank &);
    bool asend_status(const Particle_Buffer<PT> &);
    bool arecv_status(const Particle_Buffer<PT> &);
    void asend_end(const Particle_Buffer<PT> &);
    void arecv_end(const Particle_Buffer<PT> &);

    // Accessor functions.
    int num_send_nodes() const { return send_nodes.size(); }
    int num_recv_nodes() const { return recv_nodes.size(); }
    int num_bound_cells() const { return boundary_node.size(); }
    const vf_int& get_b_node() const { return boundary_node; }
    const vf_int& get_b_cell() const { return boundary_cell; }
    const sf_int& get_recv_nodes() const { return recv_nodes; }
    const sf_int& get_send_nodes() const { return send_nodes; }

    // Checks on state of buffers.
    int get_send_size() const;
    int get_recv_size() const;

    // Diagnostic functions.
    void print(std::ostream &) const;
    void arecv_print_status(const Particle_Buffer<PT> &);
};

//---------------------------------------------------------------------------//
// inline functions for the Communicator
//---------------------------------------------------------------------------//
// convert a global node index to a local node index

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

} // end namespace rtt_imc

#endif                          // __imc_Communicator_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Communicator.hh
//---------------------------------------------------------------------------//
