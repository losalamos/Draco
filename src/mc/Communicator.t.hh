//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Communicator.t.hh
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

namespace rtt_mc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for the communicator class.
 *
 * \param r_nodes vector of nodes the processor receives from
 *
 * \param s_nodes vector of nodes the processor sends to
 *
 * \param b_node  vector<vector> giving a vector of nodes for each boundary
 * cell; indexed [boundary_cell][node]
 *
 * \param b_cell  vector<vector> giving the local cell indices for a boundary
 * cell; indexed [boundary_cell][local_cell]
 */
template<class PT>
Communicator<PT>::Communicator(const sf_int &r_nodes,
			       const sf_int &s_nodes,
			       const vf_int &b_node,
			       const vf_int &b_cell)
    : recv_buffer(r_nodes.size()),
      send_buffer(s_nodes.size()), 
      recv_nodes(r_nodes), 
      send_nodes(s_nodes), 
      boundary_node(b_node), 
      boundary_cell(b_cell),
      last_node(b_cell.size())
{
    std::fill(last_node.begin(), last_node.end(), 0);
    Ensure (boundary_node.size() == boundary_cell.size());
    Ensure (recv_buffer.size() == recv_nodes.size());
    Ensure (send_buffer.size() == send_nodes.size());
}

//---------------------------------------------------------------------------//
// PARTICLE BUFFERING AND COMMUNICATION
//---------------------------------------------------------------------------//
/*!
 * \brief Buffer and send particles to other nodes.
 *
 * This function takes a particle and buffers it in an appropriate send
 * buffer for the nodes with which this processor communicates.  If the
 * buffer is full it is sent out using non-blocking communication; however,
 * we wait on the posted non-blocking sends to finish before leaving this
 * function. 
 *
 * The particle cell index \b must be negative, indicating that the particle
 * is currently in a boundary cell on this processor.
 *
 * \param particle rtt_dsxx::SP to a particle in a boundary cell
 *
 * \return the global node index if data was sent during the call; otherwise
 * return -1 if no message-passing was performed
 */
template<class PT>
int Communicator<PT>::communicate(SP_PT particle)
{
    // determine the boundary cell that this particle is in
    int bcell = -particle->get_cell();
    Require (bcell > 0);

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

    // fill the appropriate buffer
    particle->set_cell(cell);
    send_buffer[local_node].buffer_particle(*particle);

    // if buffer is full send it to the destination
    if (send_buffer[local_node].is_full()) 
    {
	// send buffer using non-blocking + wait; the buffer should
	// automatically empty itself after the wait
	send_buffer[local_node].post_asend(global_node);
	send_buffer[local_node].async_wait();
	Check (send_buffer[local_node].is_empty());
	Check (!send_buffer[local_node].comm_status());

	// return the node this particle has been sent to
	return global_node;
    }

    // if we didn't send anything out, return -1 as indicator
    return -1;
}

//---------------------------------------------------------------------------//
// MESSAGE PASSING OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Post receives to all receiving buffers that are not already
 * active. 
 */
template<class PT>
void Communicator<PT>::post()
{   
    // global node index
    int global_node = 0;

    // loop through receiving comm_buffers and post receives
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// get the global receive node index
	global_node = recv_nodes[i];
	
	// post receives if receive buffer is not active
	recv_buffer[i].post_arecv(global_node);
    }
}

//---------------------------------------------------------------------------//
// check to see if any buffers have come in, if they have stuff them in the
// bank and post a new receive for that buffer
/*!
 * \brief Check if buffers have been received and then post a new receive.
 *
 * This function checks to see if any posted receives have completed.  If so,
 * a new receive is posted.
 *
 * \return number of particles received and added to bank
 */
template<class PT> 
int Communicator<PT>::arecv_post(Bank &bank)
{
    // tag to check if something arrives from somewhere
    int arrived = 0;

    // global node index
    int global_node = 0;

    // loop over processors
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// get the global receive node index
	global_node = recv_nodes[i];
	
	// check to see if anything has come in on that node
	if (recv_buffer[i].async_check())
	{
	    Check (!recv_buffer[i].comm_status());

	    // add the number of particles to the running tally
	    arrived += recv_buffer[i].get_num_particles_in_buffer();

	    // add the data to the bank
	    recv_buffer[i].add_to_bank(bank);
	    Check (recv_buffer[i].is_empty());

	    // post a new receive to that buffer
	    recv_buffer[i].post_arecv(global_node);
	}
    }
    
    // return tag indicating where the bank has been filled
    return arrived;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait on receive buffers without posting new receives.
 *
 * Particles that are received are added to the bank.  This function waits on
 * every receive buffer, so if nothing is sent to the receive buffer this
 * will cause a hang.
 *
 * \return number of particles received and added to bank
 */
template<class PT>
int Communicator<PT>::arecv_wait(Bank &bank)
{
    // tag to check if something arrives from somewhere
    int arrived = 0;

    // loop over processors that we are receiving from
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// wait and receive the buffer
	recv_buffer[i].async_wait();

	// add the number of particles to the running tally
	arrived += recv_buffer[i].get_num_particles_in_buffer();

	// add the particles to the bank
	recv_buffer[i].add_to_bank(bank);
	Check (recv_buffer[i].is_empty());
    }
    
    // return number of particles added to the bank
    return arrived;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Receive an empty buffer without posting a new receive on all
 * receive nodes.
 *
 * This function, paired with asend_end, is used to terminate asynchronous
 * communication because no new receives are posted after the wait.  The
 * buffer that is received \b must be empty, which is the case if asend_end
 * is used.
 */
template<class PT>
void Communicator<PT>::arecv_end()
{
    // loop over processors that we are receiving from
    for (int i = 0; i < recv_nodes.size(); i++)
    {
	// wait and receive the buffer, it better have zeroes in it
	recv_buffer[i].async_wait();
	Insist (recv_buffer[i].is_empty(),
		"arecv_end received a non-empty buffer.");
	Check (!recv_buffer[i].comm_status());
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send an empty buffer to all send nodes.
 *
 * This function, paired with arecv_end, is used to terminate asynchronous
 * communication.  The send buffers \b must be empty when calling this
 * function. 
 */
template<class PT>
void Communicator<PT>::asend_end()
{
    // global node index
    int global_node = 0;

    // loop through processors and send out whatever is in there
    for (int i = 0; i < send_buffer.size(); i++)
    {
	// get the global send node index
	global_node = send_nodes[i];

	// send out asynchronously, it better have zeroes in it
	Insist (send_buffer[i].is_empty(), "asend_end buffer not empty.");
	send_buffer[i].post_asend(global_node);
	send_buffer[i].async_wait();
	Check (!send_buffer[i].comm_status());
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send out all buffers that contain data.
 *
 * Empty buffers are not sent.
 *
 * \return vector of processors that data is sent to
 */
template<class PT>
typename Communicator<PT>::sf_int Communicator<PT>::flush()
{
    // tag to check where we sent things
    sf_int sent;
    
    // global node index
    int global_node = 0;

    // loop through processors and send out whatever is in there
    for (int i = 0; i < send_buffer.size(); i++)
    {
	// get the global send node index and initialize sent to zero
	global_node = send_nodes[i];

	// send out asynchronously
	if (!send_buffer[i].is_empty())
	{
	    send_buffer[i].post_asend(global_node);
	    send_buffer[i].async_wait();
	    Check (send_buffer[i].is_empty());
	    Check (!send_buffer[i].comm_status());

	    // add the global node to the sent vector
	    sent.push_back(global_node);
	}
    }

    // return tag indicating where things have been sent
    return sent;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send out all buffers.
 *
 * All buffers are sent out whether or not they contain data.
 */
template<class PT>
void Communicator<PT>::flush_all()
{
    // global node index
    int global_node = 0;

    // loop through processors and send out whatever is in there
    for (int i = 0; i < send_buffer.size(); i++)
    {
	// get the global send node index
	global_node = send_nodes[i];

	// send out asynchronously
	send_buffer[i].post_asend(global_node);
	send_buffer[i].async_wait();
	Check (send_buffer[i].is_empty());
	Check (!send_buffer[i].comm_status());
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Free the receive buffers.
 *
 * The receive buffers \b must be empty before they are freed.
 */
template<class PT>
void Communicator<PT>::free()
{
    // loop through receive processors and free the buffers so that we can
    // move onward
    for (int i = 0; i < recv_buffer.size(); i++)
    {
	Insist (recv_buffer[i].is_empty(),
		"Tried to free a non-empty receive buffer.");

	recv_buffer[i].async_free();
	Ensure (!recv_buffer[i].comm_status());
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the status (in use or free) of the receive buffers.
 *
 * This function checks to see if \b all receive buffers are on.  If any
 * receive buffer is off then false is returned.
 */
template<class PT>
bool Communicator<PT>::arecv_status()
{
    for (int i = 0; i < recv_buffer.size(); i++)
	if (!recv_buffer[i].comm_status())
	    return false;
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the status (in use or free) of the send buffers.
 *
 * This function checks to see if \b all send buffers are on.  If any send
 * buffer is off then false is returned.
 */
template<class PT>
bool Communicator<PT>::asend_status()
{
    for (int i = 0; i < send_buffer.size(); i++)
	if (!send_buffer[i].comm_status())
	    return false;
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the status of the receive buffers to standard error.
 */
template<class PT>
void Communicator<PT>::arecv_print_status()
{
    using C4::node;
    using std::cerr;
    using std::endl;

    for (int i = 0; i < recv_buffer.size(); i++)
	cerr << "Receive buffer on node " << node() << " from node " 
	     << recv_nodes[i] << " has status "	  
	     << recv_buffer[i].comm_status() << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of particles in the send buffers.
 */
template<class PT>
int Communicator<PT>::get_send_size() const
{
    // loop through send buffers and add up particles in there
    int number = 0;
    for (int i = 0; i < send_buffer.size(); i++)
	number += send_buffer[i].get_num_particles_in_buffer();
    
    // return the total number of particles in the send buffers
    return number;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of particles in the receive buffers.
 */
template<class PT>
int Communicator<PT>::get_recv_size() const
{
    // loop through send buffers and add up particles in there
    int number = 0;
    for (int i = 0; i < recv_buffer.size(); i++)
	number += recv_buffer[i].get_num_particles_in_buffer();
    
    // return the total number of particles in the send buffers
    return number;
}

//---------------------------------------------------------------------------//
// PRINT DIAGNOSTIC FUNCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Print out the communicator for diagnostics.
 *
 * \param out ostream object to stream output to
 */
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

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Communicator.t.hh
//---------------------------------------------------------------------------//
