//----------------------------------*-C++-*----------------------------------//
// Parallel_Builder.hh
// Thomas M. Evans
// Tue Apr 14 14:50:21 1998
//---------------------------------------------------------------------------//
// @> Parallel_Builder header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Parallel_Builder_hh__
#define __imctest_Parallel_Builder_hh__

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

#include "imctest/Names.hh"
#include "imctest/Coord_sys.hh"
#include "imctest/Layout.hh"
#include "imctest/Opacity.hh"
#include "imctest/Source_Init.hh"
#include "imctest/Particle_Buffer.hh"
#include "ds++/SP.hh"

IMCSPACE

template<class MT> class Source;

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
    template<class PT>
    void dist_census(const Source_Init<MT> &, const Particle_Buffer<PT> &);

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

  // source functions
    template<class PT> 
    void send_Source(const Source_Init<MT> &, const Particle_Buffer<PT> &);
    SP<Source<MT> > recv_Source();

  // Mesh passing functionality
    void send_Mesh(const MT &);
    SP<MT> recv_Mesh();

  // Opacity passing functionality
    void send_Opacity(const Opacity<MT> &);
    SP<Opacity<MT> > recv_Opacity(SP<MT>);
};

CSPACE

#endif                          // __imctest_Parallel_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Parallel_Builder.hh
//---------------------------------------------------------------------------//
