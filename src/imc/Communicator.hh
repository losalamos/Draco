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

template<class PT>
class Communicator
{
private:
  // Comm_Buffers for particle tranport
    typename Particle_Buffer<PT>::Comm_Vector recv_buffer;
    typename Particle_Buffer<PT>::Comm_Vector send_buffer;

  // data with boundary cell->across processor info 
    vector<int> procs;
    vector<vector<int> > bcell_proc;
    vector<vector<int> > lcell_proc;

public:
  // constructors
    Communicator(int = 1);

  // buffer a particle destined for a new processor
    bool communicate(SP<PT>);

  // accessor functions
    int get_procs() const { return procs.size(); }

  // <<CONTINUE HERE>>


CSPACE

#endif                          // __imc_Communicator_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Communicator.hh
//---------------------------------------------------------------------------//
