//----------------------------------*-C++-*----------------------------------//
// Parallel_Builder.cc
// Thomas M. Evans
// Tue Apr 14 14:50:22 1998
//---------------------------------------------------------------------------//
// @> Parallel_Builder implementation file
//---------------------------------------------------------------------------//

#include "imctest/Parallel_Builder.hh"
#include "imctest/XYCoord_sys.hh"
#include "imctest/XYZCoord_sys.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <string>
#include <vector>
#include <iostream>

IMCSPACE

// C4 necessities
using C4::node;
using C4::nodes;
using C4::Send;
using C4::Recv;

// std necessities
using std::vector;
using std::string;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// Mesh passing interface
//---------------------------------------------------------------------------//
// send out the Mesh components

template<class MT>
void Parallel_Builder<MT>::send_Mesh(const MT &mesh)
{
  // send the mesh components out to the other processors

  // assure that we are on the host node only
    Check (!node());

  // let us pass a coordinate system, shall we
    send_Coord(mesh.get_Coord());

  // let us pass the Layout
    send_Layout(mesh.get_Layout());
}

//---------------------------------------------------------------------------//
// receive the Mesh components

template<class MT>
SP<MT> Parallel_Builder<MT>::recv_Mesh()
{
  // receive the Mesh components from the host and build the new mesh on this 
  // topology

  // assure that we are not on the host node
    Check (node());

  // get coordinate system
    SP<Coord_sys> coord = recv_Coord();

  // get Layout
    Layout layout(0);

  // get vertices and cell_pair
    typename MT::CCVF_a vertex(0);
    typename MT::CCVF_i cell_pair(0);

  // build mesh
    SP<OS_Mesh> mesh = new OS_Mesh(coord, layout, vertex, cell_pair);

  // return mesh
    return mesh;
}

//---------------------------------------------------------------------------//
// Mesh passing implementation
//---------------------------------------------------------------------------//
// pass the Coord_sys object

template<class MT>
void Parallel_Builder<MT>::send_Coord(const Coord_sys &coord)
{
  // send out the coordinate system designator

  // send variables
    string cs          = coord.get_Coord();
    const char *sendcs = cs.c_str();
    int cs_size        = cs.size();

  // send the Coord_sys string
    for (int np = 1; np < nodes(); np++)
    {
	Send (cs_size, np);
	Send (sendcs, cs_size, np);
    }
}

//---------------------------------------------------------------------------//
// receive the Coord_sys

template<class MT>
SP<Coord_sys> Parallel_Builder<MT>::recv_Coord()
{
  // receive and build the coordinate system

  // define Coord_sys SP
    SP<Coord_sys> coord;

  // receive the designation
    int cs_size;
    Recv (cs_size, 0);
    char *reccs = new char[cs_size];
    Recv (reccs, cs_size, 0);
    string cs = reccs;
    delete [] reccs;

  // build the coordinate system
    if (cs == "xy")
    {
	SP<XYCoord_sys> xycoord = new XYCoord_sys;
	coord = xycoord;
    }
    else if (cs == "xyz")
    {
	SP<XYZCoord_sys> xyzcoord = new XYZCoord_sys;
	coord = xyzcoord;
    }
    
  // return base class SP to a derived Coord_sys
    return coord;
}

//---------------------------------------------------------------------------//
// pass the Layout

template<class MT>
void send_Layout(const Layout &layout)
{
  // set the Layout size
    const int num_cells = layout.num_cells();

  // calculate the number of faces on each cell and total size of the Layout
  // (layout_size = SUM_(i)^(num_cells) (num_faces[i]))
    int num_faces[num_cells] = {0};
    int layout_size = 0;
    for (int i = 1; i <= num_cells; i++)
    {
	num_faces[i-1] = layout.num_faces(i);
	layout_size += num_faces[i-1];
    }

  // write the faces array for passing
    int *faces = new int[layout_size];
    int index = 0;
    for (int i = 1; i <= num_cells; i++)
	for (int j = 1; j <= layout.num_faces(i); j++)
	{
	    faces[index] = layout(i,j);
	    index++;
	}
    
  // pass all this good stuff to the other processors (for now)
    for (int np = 1; np < nodes(); np++)
    {
      // pass the Layout cell size
	Send (num_cells, np);

      // pass the Layout total size
	Send (layout_size, np);

      // pass the num_faces array
	Send (&num_faces[0], num_cells, np);

      // pass the face-values array
	Send (faces, layout_size, np);
    }

  // delete dynamically allocated faces array
    delete [] faces;
}
   

CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Builder.cc
//---------------------------------------------------------------------------//
