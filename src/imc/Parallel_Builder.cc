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

  // let us pass the vertex
    send_vertex(mesh.get_vertex());

  // let us pass the cell_pair
    send_cellpair(mesh.get_cell_pair());
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
    Layout layout = recv_Layout();

  // get vertices and cell_pair
    typename MT::CCVF_a vertex = recv_vertex();
    typename MT::CCVF_i cell_pair = recv_cellpair();

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
	Send (cs_size, np, 1);
	Send (sendcs, cs_size, np, 2);
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
    Recv (cs_size, 0, 1);
    char *reccs = new char[cs_size];
    Recv (reccs, cs_size, 0, 2);
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
void Parallel_Builder<MT>::send_Layout(const Layout &layout)
{
  // set the Layout size
    int num_cells = layout.num_cells();

  // calculate the number of faces on each cell and total size of the Layout
  // (layout_size = SUM_(i)^(num_cells) (num_faces[i]))
    int *num_faces = new int[num_cells];
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
	Send (num_cells, np, 3);

      // pass the Layout total size
	Send (layout_size, np, 4);

      // pass the num_faces array
	Send (&num_faces[0], num_cells, np, 5);

      // pass the face-values array
	Send (faces, layout_size, np, 6);
    }

  // delete dynamically allocated faces array
    delete [] faces;
    delete [] num_faces;
}

//---------------------------------------------------------------------------//
// receive the Layout

template<class MT>
Layout Parallel_Builder<MT>::recv_Layout()
{
  // receive and build the Layout

  // get the number of cells and size of the Layout
    int num_cells;
    int layout_size;
    Recv (num_cells, 0, 3);
    Recv (layout_size, 0, 4);

  // make the Layout
    Layout layout = num_cells;

  // receive the Layout data arrays
    int *num_faces = new int[num_cells];
    int *faces = new int[layout_size];
    Recv (num_faces, num_cells, 0, 5);
    Recv (faces, layout_size, 0, 6);

  // rebuild the Layout
    int index = 0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	layout.set_size(cell, num_faces[cell-1]);
	for (int i = 1; i <= num_faces[cell-1]; i++)
	{
	    layout(cell, i) = faces[index];
	    index++;
	}
    }

  // get rid of dynamic arrays
    delete [] num_faces;
    delete [] faces;
    
  // return the new Layout on this node
    return layout;
}

//---------------------------------------------------------------------------//
// pass the mesh vertices

template<class MT>
void Parallel_Builder<MT>::send_vertex(const typename MT::CCVF_a &vertex)
{
  // send the vertices to another processor
    
  // get all necessary vertex dimensions
    int vertex_dim  = vertex.size();
    int vertex_size = vertex[0].size();
    int total_size  = vertex_dim * vertex_size;

  // push all vertices onto one array for communication
    double *vert = new double[total_size];
    int index = 0;
    for (int d = 0; d < vertex_dim; d++)
    {
      // some assertions to check that the vertex size is constant over all
      // dimensions 
	Check (vertex_size == vertex[d].size());

      // now load up the communication vertex array
	for (int i = 0; i < vertex_size; i++)
	{
	    vert[index] = vertex[d][i];
	    index++;
	}
    }

  // send out the goodies
    for (int np = 1; np < nodes(); np++)
    {
      // send the dimensionality of the vertex-array
	Send (vertex_dim, np, 7);

      // send the total size of the vertex-array
	Send (total_size, np, 8);

      // send the vertex array
	Send (vert, total_size, np, 9);
    }

  // delete dynamically allocated arrays
    delete [] vert;
}
	   
//---------------------------------------------------------------------------//
// receive the vertex

template<class MT>
typename MT::CCVF_a Parallel_Builder<MT>::recv_vertex()
{
  // receive and rebuild the vertices

  // first get the sizes
    int vertex_dim;
    int total_size;
    Recv (vertex_dim, 0, 7);
    Recv (total_size, 0, 8);

  // find the size of each dimensional vertex array
    Check (!(total_size % vertex_dim));
    int vertex_size = total_size / vertex_dim;

  // make a new CCVF_a vertex array
    typename MT::CCVF_a vertex(vertex_dim);

  // get the vertices from the host
    double *vert = new double[total_size];
    Recv (vert, total_size, 0, 9);

  // assign values to the new vertex object
    int index = 0;
    for (int d = 0; d < vertex_dim; d++)
    {
	vertex[d].resize(vertex_size);
	for (int i = 0; i < vertex_size; i++)
	{
	    vertex[d][i] = vert[index];
	    index++;
	}
    }

  // delete dynamic allocated arrays
    delete [] vert;

  // return the new vertex guy
    return vertex;
}

//---------------------------------------------------------------------------//
// pass the cell_pair

template<class MT>
void Parallel_Builder<MT>::send_cellpair(const typename MT::CCVF_i
					 &cell_pair)
{
  // send the cell_pair array to another processor
    
  // set the cell_pair size
    int num_cells = cell_pair.size();

  // calculate the number of vertices per cell and the total size of the
  // cell_pair object
    int *num_vert = new int[num_cells];
    int size = 0;
    for (int i = 0; i < num_cells; i++)
    {
	num_vert[i] = cell_pair[i].size();
	size += num_vert[i];
    }

  // write the cell_pair array for passing
    int *vertices = new int[size];
    int index = 0;
    for (int i = 0; i < num_cells; i++)
	for (int j = 0; j < cell_pair[i].size(); j++)
	{
	    vertices[index] = cell_pair[i][j];
	    index++;
	}

  // pass the important stuff to other processors
    for (int np = 1; np < nodes(); np++)
    {
      // pass the size
	Send (num_cells, np, 10);

      // pass the total size of the cell_pair object
	Send (size, np, 11);

      // pass the num_vert array
	Send (num_vert, num_cells, np, 12);
	
      // pass the vertices-values array
	Send (vertices, size, np, 13);
    }

  // delete the dynamically allocated arrays
    delete [] num_vert;
    delete [] vertices;
}

//---------------------------------------------------------------------------//
// receive the cell_pair object

template<class MT>
typename MT::CCVF_i Parallel_Builder<MT>::recv_cellpair()
{
  // receive and rebuild the cell_pair array
    
  // first get the sizes 
    int num_cells;
    int size;
    Recv (num_cells, 0, 10);
    Recv (size, 0, 11);

  // make a new cell_pair object
    typename MT::CCVF_i cell_pair(num_cells);

  // receive the cell_pair data
    int *num_vert = new int[num_cells];
    int *vertices = new int[size];
    Recv (num_vert, num_cells, 0, 12);
    Recv (vertices, size, 0, 13);

  // rebuild the cell_pair object
    int index = 0;
    for (int i = 0; i < num_cells; i++)
    {
	cell_pair[i].resize(num_vert[i]);
	for (int j = 0; j < num_vert[i]; j++)
	{
	    cell_pair[i][j] = vertices[index];
	    index++;
	}
    }
    Check (index == size);

  // delete dynamic arrays
    delete [] num_vert;
    delete [] vertices;

  // return the cell_pair object
    return cell_pair;
}
    
CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Builder.cc
//---------------------------------------------------------------------------//
