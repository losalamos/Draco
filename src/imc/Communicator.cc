//----------------------------------*-C++-*----------------------------------//
// Communicator.cc
// Thomas M. Evans
// Thu Jul  9 19:16:28 1998
//---------------------------------------------------------------------------//
// @> Communicator class implementation file
//---------------------------------------------------------------------------//

#include "imc/Communicator.hh"

IMCSPACE

// std components
using std::fill;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// general constructor

template<class PT>
Communicator<PT>::Communicator(const vector<int> &r_nodes,
			       const vector<int> &s_nodes,
			       const vector<vector<int> > &b_node,
			       const vector<vector<int> > &b_cell)
    : recv_buffer(r_nodes.size()), send_buffer(s_nodes.size()), 
      recv_nodes(r_nodes), send_nodes(s_nodes), boundary_node(b_node), 
      boundary_cell(b_cell), last_node(b_cell.size())
{
    fill(last_node.begin(), last_node.end(), 0);
    Ensure (boundary_node.size() == boundary_cell.size());
}

//---------------------------------------------------------------------------//
// particle buffering and communication
//---------------------------------------------------------------------------//
// take a particle and place it in the appropriate buffer to send to another
// node

template<class PT>
int Communicator<PT>::communicate(const Particle_Buffer<PT> &buffer,
				  SP<PT> particle)
{
  // determine the boundary cell that this particle is in
    int bcell = -particle->get_cell();

  // determine the index for the destination processor and the local cell 
  // index of bcell on that processor
    last_node[bcell-1] = (last_node[bcell-1] < boundary_node[bcell-1].size()) 
	? last_node[bcell-1] : 0;
    int index = last_node[bcell-1]++;

  // get the local node id and the local cell id and the global node id
    int node        = boundary_node[bcell-1][index];
    int cell        = boundary_cell[bcell-1][index];
    int global_node = send_nodes[node];

  // fill the appropriate buffer
    particle->set_cell(cell);
    buffer.buffer_particle(send_buffer[node], *particle);

  // if buffer is full send it to the destination
    if (send_buffer[node].n_part == Particle_Buffer<PT>::get_buffer_s())
    {
	buffer.asend_buffer(send_buffer[node], global_node);
	Check (send_buffer[node].n_part == 0);
    }

  // return the node this particle has been sent to
    return global_node;
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
// bank

template<class PT>
bool Communicator<PT>::arecv(const Particle_Buffer<PT> &buffer, 
			     typename Particle_Buffer<PT>::Bank &bank)
{
  // tag to check if something arrives from somewhere
    int arrived = 0;

  // loop over processors
    for (int i = 0; i < recv_nodes.size(); i++)
    {
      // get the global receive node index
	int global_node = recv_nodes[i];

      // check to see if anything has come in on that node
	if (buffer.check_arecv(recv_buffer[i], global_node))
	{
	    buffer.add_to_bank(recv_buffer[i], bank);
	    Check (recv_buffer[i].n_part == 0);
	    arrived++;
	}
    }
    
  // return tag indicating if the bank has been filled
    return arrived;
}

//---------------------------------------------------------------------------//
// flush out the send buffers

template<class PT>
void Communicator<PT>::flush(const Particle_Buffer<PT> &buffer)
{
  // loop through processors and send out whatever is in there
    for (int i = 0; i < send_buffer.size(); i++)
    {
      // get the global send node index
	int global_node = send_nodes[i];

      // send out asynchronously
	buffer.asend_buffer(send_buffer[i], global_node);
	Check (send_buffer[i].n_part == 0);
    }
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

CSPACE

//---------------------------------------------------------------------------//
//                              end of Communicator.cc
//---------------------------------------------------------------------------//
