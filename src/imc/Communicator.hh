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

IMCSPACE

// DRACO components
using dsxx::SP;

// std components
using std::vector;

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

  // checks on state of buffers
    int get_send_size() const;
    int get_recv_size() const;

  // need a function to take a global_node index and convert it to a local
  // node index!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
};

//---------------------------------------------------------------------------//
// inline functions for the Communicator
//---------------------------------------------------------------------------//

CSPACE

#endif                          // __imc_Communicator_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Communicator.hh
//---------------------------------------------------------------------------//
