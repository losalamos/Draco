//----------------------------------*-C++-*----------------------------------//
// Parallel_Builder.hh
// Thomas M. Evans
// Tue Apr 14 14:50:21 1998
//---------------------------------------------------------------------------//
// @> Parallel_Builder header file
//---------------------------------------------------------------------------//

#ifndef __imc_Parallel_Builder_hh__
#define __imc_Parallel_Builder_hh__

//===========================================================================//
// class Parallel_Builder - 
//
// Purpose : pass the mesh, opacity, and source info to IMC_topology
//           processors
//
// revision history:
// -----------------
//  0) original
//  1)  7-28-98 : fixed calculation of cell_pair data in build_cells() member 
//                function; this error did not cause any transport errors 
//                luckily
//  2)  9-16-98 : added global_cells vector which maps the local cells to
//                global cells, this data is sent to IMC processors in the 
//                constructors for Parallel_Builder; added some DBC Requires
//                to mesh topology query functions
//  3)  9-16-98 : in dist_census, added global_cell -> local_cell update for
//                particles that are sent out
//  4)  10-9-98 : added submesh=true indicator to send and recv_Mesh
//                functions; Meshes that are sent out should be designated 
//                submeshes.
//  5) 10-20-98 : fixed array index in recv_Communicator; changed
//                recv_nodes[node_size + 1] to recv_nodes[node_size + i]
// 
//===========================================================================//

#include "mc/Coord_sys.hh"
#include "Layout.hh"
#include "Opacity.hh"
#include "Source_Init.hh"
#include "Particle_Buffer.hh"
#include "Source.hh"
#include "Mat_State.hh"
#include "Communicator.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

namespace rtt_imc 
{

// draco stuff
using C4::node;
using C4::nodes;
using rtt_rng::Rnd_Control;
using dsxx::SP;

// std stuff
using std::vector;
using std::find;
using std::string;

template<class MT>
class Parallel_Builder
{
private:
  // topology data
    vector<vector<int> > cells_per_proc;
    vector<vector<int> > procs_per_cell;
    vector<vector<int> > bound_cells;
    vector<int> global_cells;

  // decomposition mode
    string parallel_scheme;

  // calculate topology map
    void parallel_topology(const MT &, const Source_Init<MT> &);

  // calculate parallel source distributions

  // distribute the census, volume, and surface sources
    template<class PT>
    void dist_census(const Source_Init<MT> &, const Particle_Buffer<PT> &,
		     typename Particle_Buffer<PT>::Census &);
    void dist_vol(const Source_Init<MT> &, typename MT::CCSF_int &,
		  typename MT::CCSF_int &, typename MT::CCSF_double &,
		  typename MT::CCVF_double &);
    void dist_ss(const Source_Init<MT> &, typename MT::CCSF_int &,
		 typename MT::CCSF_int &, typename MT::CCSF_int &,
		 typename MT::CCSF_double &);

  // receive the census, volume, and surface sources
    template<class PT>
    void recv_census(const Particle_Buffer<PT> &,
		     typename Particle_Buffer<PT>::Census &);
    void recv_vol(typename MT::CCSF_int &, typename MT::CCSF_int &,
		  typename MT::CCSF_double &, typename MT::CCVF_double &);
    void recv_ss(typename MT::CCSF_int &, typename MT::CCSF_int &,
		 typename MT::CCSF_int &, typename MT::CCSF_double &);

  // functionality for Mesh passing

  // build and pass the Layout
    void send_Layout(rtt_mc::Layout &);
    rtt_mc::Layout recv_Layout();
    rtt_mc::Layout build_Layout(const rtt_mc::Layout &, int);

  // pass the coord_sys
    void send_Coord(SP<rtt_mc::Coord_sys>);
    SP<rtt_mc::Coord_sys> recv_Coord();
  
  // pass the cells (vertices and cell pairings)
  // overload this function to allow usage with both regular OS meshes and
  // the CAR meshes (the latter needs to send the generation also)
    void send_cells(const MT &, typename MT::NCVF_d &, typename MT::CCVF_i &);
    void send_cells(const MT &, typename MT::NCVF_d &, typename MT::CCVF_i &, 
		    typename MT::CCSF_i &);
  // overload this function to allow usage with both regular OS meshes and
  // the CAR meshes (the latter needs to receive the generation also)
    void build_cells(const MT &, typename MT::NCVF_d &, typename MT::CCVF_i &,
		     int);
    void build_cells(const MT &, typename MT::NCVF_d &, typename MT::CCVF_i &,
		     typename MT::CCSF_i &, int);
    void send_vertex(const typename MT::NCVF_d &, int);
    typename MT::NCVF_d recv_vertex();
    void send_cellpair(const typename MT::CCVF_i &, int);
    typename MT::CCVF_i recv_cellpair();
    void send_generation(const typename MT::CCSF_i &, int);
    typename MT::CCSF_i recv_generation();

  // Communicator functionality
    template<class PT> SP<Communicator<PT> > build_Communicator(int);

public:
  // default constructor, for non-host nodes
    Parallel_Builder();
    
  // constructor for host-node
    Parallel_Builder(const MT &, const Source_Init<MT> &);

  // Mesh passing functionality
    SP<MT> send_Mesh(const MT &);
    SP<MT> recv_Mesh();

  // Opacity passing functionality
    SP<Opacity<MT> > send_Opacity(SP<MT>, const Opacity<MT> &);
    SP<Opacity<MT> > recv_Opacity(SP<MT>);

  // Mat State passing functionality
    SP<Mat_State<MT> > send_Mat(SP<MT>, const Mat_State<MT> &);
    SP<Mat_State<MT> > recv_Mat(SP<MT>);

  // source passing functionality
    template<class PT> SP<Source<MT> > 
    send_Source(SP<MT>, SP<Mat_State<MT> >, SP<Rnd_Control>, 
		const Source_Init<MT> &, const Particle_Buffer<PT> &);  
    template<class PT> SP<Source<MT> > 
    recv_Source(SP<MT>, SP<Mat_State<MT> >, SP<Rnd_Control>,
		const Particle_Buffer<PT> &); 

  // communicator mapping functionality
    template<class PT> 
    SP<Communicator<PT> > send_Communicator();
    template<class PT>
    SP<Communicator<PT> > recv_Communicator();

  // Mesh mapping functionality functions available on the host node
    inline int master_cell(int icell, int proc) const;
    inline int imc_cell(int mcell, int proc) const;
    inline vector<int> get_cells(int proc) const;
    inline vector<int> get_procs(int mcell) const;
    int num_cells() const { return procs_per_cell.size(); }
    int num_cells(int proc) const { return cells_per_proc[proc].size(); }
    int num_procs(int mcell) const { return procs_per_cell[mcell-1].size(); } 
    
  // Mesh mapping functionality available on all nodes
    int master_cell(int icell) const { return global_cells[icell-1]; }

  // parallel scheme
    string get_parallel_scheme() const { return parallel_scheme; }

  // diagnostics
    void print(ostream &) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
ostream& operator<<(ostream &out, const Parallel_Builder<MT> &object)
{
    object.print(out);
    return out;
}

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR PARALLEL_BUILDER
//---------------------------------------------------------------------------//
// determine a cell on the master mesh given a cell on a particular processor

template<class MT>
inline int Parallel_Builder<MT>::master_cell(int icell, int proc) const
{
    Require (!node());
    Require (proc < nodes());
    Require (icell <= cells_per_proc[proc].size());
    return cells_per_proc[proc][icell-1];
}

//---------------------------------------------------------------------------//
// determine a cell on the IMC mesh on processor given the master cell,
// returns 0 (false) if the cell isn't on this processor

template<class MT>
inline int Parallel_Builder<MT>::imc_cell(int mcell, int proc) const
{
    Require (!node());
    Require (proc < nodes());

  // get the iterator location, for const_iterator explanation see KAI
  // response of 6-19-98
    vector<int>::const_iterator itr = find(cells_per_proc[proc].begin(), 
					   cells_per_proc[proc].end(),
					   mcell);   

  // if the cell is here return it
    if (itr != cells_per_proc[proc].end())
	return (itr - cells_per_proc[proc].begin()) + 1;
    return 0;
}

//---------------------------------------------------------------------------//
// return a vector<int> holding the master cells on this processor

template<class MT>
inline vector<int> Parallel_Builder<MT>::get_cells(int proc) const
{
    Require (!node());
    Require (proc < nodes());
    return cells_per_proc[proc];
}

//---------------------------------------------------------------------------//
// return a vector<int> holding the procs that contain this master cell

template<class MT>
inline vector<int> Parallel_Builder<MT>::get_procs(int mcell) const
{
    Require (!node());
    return procs_per_cell[mcell-1];
}

} // end namespace rtt_imc

#endif                          // __imc_Parallel_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Parallel_Builder.hh
//---------------------------------------------------------------------------//
