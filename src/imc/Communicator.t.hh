//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Communicator.t.hh
 * \author Thomas M. Evans
 * \date   Thu Jul  9 19:16:28 1998
 * \brief  Communicator class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Communicator.hh"
#include "c4/global.hh"
#include <iomanip>

namespace rtt_imc 
{

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// general constructor

template<class PT>
Communicator<PT>::Communicator(const sf_int &r_nodes,
			       const sf_int &s_nodes,
			       const vf_int &b_node,
			       const vf_int &b_cell)
    : recv_buffer(r_nodes.size()),
      send_buffer(s_nodes.size()), 
      recv_nodes(r_nodes), send_nodes(s_nodes), 
      boundary_node(b_node), 
      boundary_cell(b_cell), last_node(b_cell.size())
{
    std::fill(last_node.begin(), last_node.end(), 0);
    Ensure (boundary_node.size() == boundary_cell.size());
    Ensure (recv_buffer.size() == recv_nodes.size());
    Ensure (send_buffer.size() == send_nodes.size());
}

//---------------------------------------------------------------------------//
// particle buffering and communication
//---------------------------------------------------------------------------//
// take a particle and place it in the appropriate buffer to send to another
// node

template<class PT>
int Communicator<PT>::communicate(const Particle_Buffer<PT> &buffer,
				  SP_PT particle)
{
    // determine the boundary cell that this particle is in
    int bcell = -particle->get_cell();
    Check (bcell > 0);

    // determine the index for the destination processor and the local cell 
    // index of bcell on that processor
    last_node[bcell-1] = (last_node[bcell-1] < boundary_node[bcell-1].size()) 
	? last_node[bcell-1] : 0;
    int index = last_node[bcell-1]++;

    // get the local node id and the local cell id and the global node id

    // convert global-to-local
    int global_node = boundary_node[bcell-1][index];
    int cell        = boundary_cell[bcell-1][index];
    int local_node  = global_to_local(global_node);
    int multiplier  = -1;

    // fill the appropriate buffer
    particle->set_cell(cell);
    buffer.buffer_particle(send_buffer[local_node], *particle);

    // if buffer is full send it to the destination
    if (send_buffer[local_node].n_part ==
	Particle_Buffer<PT>::get_buffer_s()) 
    {
	buffer.asend_buffer(send_buffer[local_node], global_node);
	buffer.async_wait(send_buffer[local_node]);
	Check (send_buffer[local_node].n_part == 0);

	// return the node this particle has been sent to
	return global_node;
    }

    // if we didn't send anything out, return -1 as indicator
    return multiplier;
}

//---------------------------------------------------------------------------//
// message passing operations
//---------------------------------------------------------------------------//
// post async receives from all recv'ing nodes

template<class PT>
void Communicator<PT>::post(const Particle_Buffer<PT> &buffer)
{
    // loop through receiving comm_buffers and post receives
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// get the global receive node index
	int global_node = recv_nodes[i];
	
	// post receives
	buffer.post_arecv(recv_buffer[i], global_node);
    }
}

//---------------------------------------------------------------------------//
// check to see if any buffers have come in, if they have stuff them in the
// bank and post a new receive for that buffer

template<class PT> 
int Communicator<PT>::arecv_post(const Particle_Buffer<PT> &buffer, 
				 Bank &bank)
{
    // tag to check if something arrives from somewhere
    int arrived = 0;

    // loop over processors
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// get the global receive node index
	int global_node = recv_nodes[i];

	// check to see if anything has come in on that node
	if (buffer.async_check(recv_buffer[i]))
	{
	    // receive the data
	    buffer.async_wait(recv_buffer[i]);
	    buffer.add_to_bank(recv_buffer[i], bank);
	    Check (recv_buffer[i].n_part == 0);
	    arrived++;

	    // post a new receive to that buffer
	    buffer.post_arecv(recv_buffer[i], global_node);
	}
    }
    
    // return tag indicating where the bank has been filled
    return arrived;
}

//---------------------------------------------------------------------------//
// wait on each buffer to receive something, stuff them in a bank, DO NOT
// post a new receive

template<class PT>
void Communicator<PT>::arecv_wait(const Particle_Buffer<PT> &buffer, 
				  Bank &bank)
{
    // loop over processors that we are receiving from
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// get the global receive node index
	int global_node = recv_nodes[i];

	// wait and receive the buffer
	buffer.async_wait(recv_buffer[i]);
	buffer.add_to_bank(recv_buffer[i], bank);
	Check (recv_buffer[i].n_part == 0);
    }
}

//---------------------------------------------------------------------------//
// wait on each buffer to receive a Comm_Buffer of zeroes, this is used to
// terminate asynchronous communication, it does not post a new receive, it
// is meant to work with asend_end

template<class PT>
void Communicator<PT>::arecv_end(const Particle_Buffer<PT> &buffer)
{
    // loop over processors that we are receiving from
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// get the global receive node index
	int global_node = recv_nodes[i];

	// wait and receive the buffer, it better have zeroes in it
	buffer.async_wait(recv_buffer[i]);
	Check (recv_buffer[i].n_part == 0);
    }
}

//---------------------------------------------------------------------------//
// send out all buffers, they should be full of zeroes, this is used to
// terminate asynchronous commuication, it is meant to work with arecv_end

template<class PT>
void Communicator<PT>::asend_end(const Particle_Buffer<PT> &buffer)
{
    // loop through processors and send out whatever is in there
    for (int i = 0; i < send_buffer.size(); i++)
    {
	// get the global send node index
	int global_node = send_nodes[i];

	// send out asynchronously, it better have zeroes in it
	Check (send_buffer[i].n_part == 0);
	buffer.asend_buffer(send_buffer[i], global_node);
	buffer.async_wait(send_buffer[i]);
    }
}

//---------------------------------------------------------------------------//
// flush out the send buffers if they contain data

template<class PT>
typename Communicator<PT>::sf_int 
Communicator<PT>::flush(const Particle_Buffer<PT> &buffer)
{
    // tag to check where we sent things
    sf_int sent(send_buffer.size());

    // loop through processors and send out whatever is in there
    for (int i = 0; i < send_buffer.size(); i++)
    {
	// get the global send node index and initialize sent to zero
	int global_node = send_nodes[i];
	sent[i] = -1;

	// send out asynchronously
	if (send_buffer[i].n_part > 0)
	{
	    buffer.asend_buffer(send_buffer[i], global_node);
	    buffer.async_wait(send_buffer[i]);
	    Check (send_buffer[i].n_part == 0);
	    sent[i] = global_node;
	}
    }

    // return tag indicating where things have been sent
    return sent;
}

//---------------------------------------------------------------------------//
// flush out all the send buffers whether they have data or are zeroes

template<class PT>
void Communicator<PT>::flush_all(const Particle_Buffer<PT> &buffer)
{
    // loop through processors and send out whatever is in there
    for (int i = 0; i < send_buffer.size(); i++)
    {
	// get the global send node index
	int global_node = send_nodes[i];

	// send out asynchronously
	buffer.asend_buffer(send_buffer[i], global_node);
	buffer.async_wait(send_buffer[i]);
	Check (send_buffer[i].n_part == 0);
    }
}

//---------------------------------------------------------------------------//
// free the receive buffers

template<class PT>
void Communicator<PT>::free(const Particle_Buffer<PT> &buffer)
{
    // loop through receive processors and free the buffers so that we can move 
    // onward
    for (int i = 0; i < recv_buffer.size(); i++)
    {
	Check (recv_buffer[i].n_part == 0);
	buffer.async_free(recv_buffer[i]);
    }
}

//---------------------------------------------------------------------------//
// get the status (in use or free) of async recv Comm_Buffers
// NOTE: this only properly checks for true conditions, that is ALL of the
// recv_buffers are active

template<class PT>
bool Communicator<PT>::arecv_status(const Particle_Buffer<PT> &buffer)
{
    for (int i = 0; i < recv_buffer.size(); i++)
	if (!buffer.comm_status(recv_buffer[i]))
	    return false;
    return true;
}

//---------------------------------------------------------------------------//
// print out the status of the recv nodes

template<class PT>
void Communicator<PT>::arecv_print_status(const Particle_Buffer<PT> &buffer)
{
    using C4::node;
    using std::cerr;
    using std::endl;

    for (int i = 0; i < recv_buffer.size(); i++)
	cerr << "Receive buffer on node " << node() << " from node " 
	     << recv_nodes[i] << " has status "	     
	     << buffer.comm_status(recv_buffer[i]) << endl;
}

//---------------------------------------------------------------------------//
// get the status (in use or free) of async send Comm_Buffers
// NOTE: this only properly checks for true conditions, that is ALL of the
// send_buffers are active

template<class PT>
bool Communicator<PT>::asend_status(const Particle_Buffer<PT> &buffer)
{
    for (int i = 0; i < send_buffer.size(); i++)
	if (!buffer.comm_status(send_buffer[i]))
	    return false;
    return true;
}

//---------------------------------------------------------------------------//
// check the number of particles in the send buffers

template<class PT>
int Communicator<PT>::get_send_size() const
{
    // loop through send buffers and add up particles in there
    int number = 0;
    for (int i = 0; i < send_buffer.size(); i++)
	number += send_buffer[i].n_part;
    
    // return the total number of particles in the send buffers
    return number;
}

//---------------------------------------------------------------------------//
// check the number of particles in the recv buffers

template<class PT>
int Communicator<PT>::get_recv_size() const
{
    // loop through send buffers and add up particles in there
    int number = 0;
    for (int i = 0; i < recv_buffer.size(); i++)
	number += recv_buffer[i].n_part;
    
    // return the total number of particles in the send buffers
    return number;
}

//---------------------------------------------------------------------------//
// print diagnostic function
//---------------------------------------------------------------------------//
// print out the communicator

template<class PT>
void Communicator<PT>::print(std::ostream &out) const
{
    using C4::node;
    using std::endl;
    using std::ios;
    using std::setw;

    out << endl;
    out << ">>> COMMUNICATOR <<<" << endl;
    out << "====================" << endl;

    out.setf(ios::right, ios::adjustfield);
    out << endl;
    out << setw(20) << "Node:" << setw(8) << node() << endl;
    out << setw(20) << "Boundary cells:" << setw(8) << boundary_node.size()
	<< endl;
    out << setw(20) << "Recv buffers:" << setw(8) << recv_nodes.size()
	<< endl;
    out << setw(20) << "Send buffers:" << setw(8) << send_nodes.size()
	<< endl;

    out << endl;
    out << setw(15) << "Boundary Cell" << setw(15) << "Send Node" << setw(15) 
	<< "Global Node" << setw(15) << "Local Cell" << endl;
    out << "------------------------------------------------------------" 
	<< endl;
    for (int i = 0; i < boundary_cell.size(); i++)
	for (int j = 0; j < boundary_cell[i].size(); j++)
	    out << setw(15) << i+1 << setw(15) << 
		global_to_local(boundary_node[i][j]) << setw(15) << 
		boundary_node[i][j] <<setw(15) << boundary_cell[i][j] <<
		endl; 

    out << endl;
    out << setw(15) << "Receive Nodes" << endl;
    out << "---------------" << endl;
    for (int i = 0; i < recv_nodes.size(); i++)
	out << setw(15) << recv_nodes[i] << endl;

    out << endl;
    out << setw(15) << "Send Nodes" << endl;
    out << "---------------" << endl;
    for (int i = 0; i < send_nodes.size(); i++)
	out << setw(15) << send_nodes[i] << endl;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Communicator.t.hh
//---------------------------------------------------------------------------//
