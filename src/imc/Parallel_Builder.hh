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
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/Coord_sys.hh"
#include "imc/Layout.hh"
#include "imc/Opacity.hh"
#include "imc/Source_Init.hh"
#include "imc/Particle_Buffer.hh"
#include "imc/Source.hh"
#include "imc/Mat_State.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>

IMCSPACE

// draco stuff
using C4::node;
using C4::nodes;
using RNG::Rnd_Control;

// std stuff
using std::vector;
using std::find;

template<class MT>
class Parallel_Builder
{
private:
  // topology data
    vector<vector<int> > cells_per_proc;
    vector<vector<int> > procs_per_cell;

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
    template<class PT> void recv_census(const Particle_Buffer<PT> &);
    void recv_vol(typename MT::CCSF_int &, typename MT::CCSF_int &,
		  typename MT::CCSF_double &, typename MT::CCVF_double &);
    void recv_ss(typename MT::CCSF_int &, typename MT::CCSF_int &,
		 typename MT::CCSF_int &, typename MT::CCSF_double &);

  // functionality for Mesh passing

  // pass the Layout
    void send_Layout(const Layout &);
    Layout recv_Layout();

  // pass the coord_sys
    void send_Coord(const Coord_sys &);
    SP<Coord_sys> recv_Coord();
  
  // pass the vertices
    void send_vertex(const typename MT::CCVF_a &);
    typename MT::CCVF_a recv_vertex();
    void send_cellpair(const typename MT::CCVF_i &);
    typename MT::CCVF_i recv_cellpair();

public:
  // default constructor, for non-host nodes
    Parallel_Builder();
    
  // constructor for host-node
    Parallel_Builder(const MT &, const Source_Init<MT> &);

  // Mesh passing functionality
    void send_Mesh(const MT &);
    SP<MT> recv_Mesh();

  // Opacity passing functionality
    void send_Opacity(const Opacity<MT> &);
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

  // Mesh mapping functionality
    inline int master_cell(int icell, int proc) const;
    inline int imc_cell(int mcell, int proc) const;
    inline vector<int> get_cells(int proc) const;
    inline vector<int> get_procs(int mcell) const;
    int num_cells() const { return procs_per_cell.size(); }
    int num_cells(int proc) const { return cells_per_proc[proc].size(); }
    int num_procs(int mcell) const { return procs_per_cell[mcell-1].size(); } 
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR PARALLEL_BUILDER
//---------------------------------------------------------------------------//
// determine a cell on the master mesh given a cell on a particular processor

template<class MT>
inline int Parallel_Builder<MT>::master_cell(int icell, int proc) const
{
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
    Require (proc < nodes());
    return cells_per_proc[proc];
}

//---------------------------------------------------------------------------//
// return a vector<int> holding the procs that contain this master cell

template<class MT>
inline vector<int> Parallel_Builder<MT>::get_procs(int mcell) const
{
    return procs_per_cell[mcell-1];
}

CSPACE

#endif                          // __imc_Parallel_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Parallel_Builder.hh
//---------------------------------------------------------------------------//
