//----------------------------------*-C++-*----------------------------------//
// Communicator.hh
// Thomas M. Evans
// Thu Jul  9 19:16:28 1998
//---------------------------------------------------------------------------//
// @> Communicator class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Communicator_hh__
#define __imc_Communicator_hh__

//===========================================================================//
// class Communicator - 
//
// Purpose : interface class which determine where a particle needs to go 
//           when it crosses a processor-domain boundary
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/Particle_Buffer.hh"
#include "ds++/SP.hh"
#include <vector>
#include <algorithm>

IMCSPACE

// DRACO components
using dsxx::SP;

// std components
using std::vector;
using std::find;

template<class PT>
class Communicator
{
private:
  // Comm_Buffers for particle tranport
    typename Particle_Buffer<PT>::Comm_Vector recv_buffer;
    typename Particle_Buffer<PT>::Comm_Vector send_buffer;

  // data with boundary cell->across processor info 
    vector<int> recv_nodes;
    vector<int> send_nodes;
    vector<int> last_node;
    vector<vector<int> > boundary_node;
    vector<vector<int> > boundary_cell;

  // convert global->local nodes
    inline int global_to_local(int) const;

  // class specific type-defs
    typedef typename Particle_Buffer<PT>::Bank Bank;

public:
  // constructors
    explicit Communicator(const vector<int> &, const vector<int> &,
			  const vector<vector<int> > &,
			  const vector<vector<int> > &);

  // buffer a particle destined for a new processor
    int communicate(const Particle_Buffer<PT> &, SP<PT>);

  // message passing (async and sync) operations
    void post(const Particle_Buffer<PT> &);
    bool arecv(const Particle_Buffer<PT> &, Bank &);
    void flush(const Particle_Buffer<PT> &);

  // accessor functions
    int num_send_nodes() const { return send_nodes.size(); }
    int num_recv_nodes() const { return recv_nodes.size(); }
    int num_bound_cells() const { return boundary_node.size(); }
    const vector<vector<int> >& get_b_node() const { return boundary_node; }
    const vector<vector<int> >& get_b_cell() const { return boundary_cell; }
    const vector<int>& get_recv_nodes() const { return recv_nodes; }
    const vector<int>& get_send_nodes() const { return send_nodes; }

  // checks on state of buffers
    int get_send_size() const;
    int get_recv_size() const;
};

//---------------------------------------------------------------------------//
// inline functions for the Communicator
//---------------------------------------------------------------------------//
// convert a global node index to a local node index

template<class PT>
int Communicator<PT>::global_to_local(int global_index) const
{
    vector<int>::const_iterator itr = find(send_nodes.begin(), 
					   send_nodes.end(), global_index);
    Ensure (itr != send_nodes.end());
    Ensure (*itr == global_index);
    return itr - send_nodes.begin();
}

CSPACE

#endif                          // __imc_Communicator_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Communicator.hh
//---------------------------------------------------------------------------//
