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
#include "ds++/SP.hh"

IMCSPACE

template<class MT>
class Parallel_Builder
{
private:

  // functionality for Mesh passing

  // pass the Layout
    void send_Layout(const Layout &);
    SP<Layout> recv_Layout();

  // pass the coord_sys
    void send_Coord(const Coord_sys &);
    SP<Coord_sys> recv_Coord();
  
  // pass the vertices
    void send_vertex(typename MT::CCVF_a &);
    SP<typename MT::CCVF_a> recv_vertex();
    void send_cellpair(typename MT::CCVF_i &);
    SP<typename MT::CCVF_a> recv_cellpair();

public:
  // default constructor
    Parallel_Builder() {}

  // Mesh passing functionality
    void send_Mesh(const MT &);
    SP<MT> recv_Mesh();
};

CSPACE

#endif                          // __imctest_Parallel_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Parallel_Builder.hh
//---------------------------------------------------------------------------//
