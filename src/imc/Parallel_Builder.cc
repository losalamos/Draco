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
#include "imctest/Math.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>

IMCSPACE

// draco necessities
using C4::node;
using C4::nodes;
using C4::Send;
using C4::Recv;
using Global::max;
using Global::min;

// std necessities
using std::vector;
using std::string;
using std::fill;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// host node constructor that determines the toplogy on the IMC-nodes

// template<class MT>
// Parallel_Builder<MT>::Parallel_Builder(const MT &mesh, 
// 				       const Source_Init<MT> &sinit)
//     : cells_per_proc(nodes()), procs_per_cell(mesh.num_cells())
// {
//   // calculate the parameters for splitting the problem amongst many
//   // processors 
//   //  parallel_params(sinit);
// }

//---------------------------------------------------------------------------//
// default constructor for IMC-nodes

template<class MT>
Parallel_Builder<MT>::Parallel_Builder()
{
  // the IMC-nodes don't use this data
}

//---------------------------------------------------------------------------//
// IMC Topology calculation
//---------------------------------------------------------------------------//
// calculate topology parameters to determine what cells go where

template<class MT>
void Parallel_Builder<MT>::parallel_params(const MT &mesh,
					   const Source_Init<MT> &sinit)
{
  // make sure we have a Source_Init and Mesh
    assert (&sinit);
    assert (&mesh);

  // number of cells in the mesh
    int num_cells = mesh.num_cells();

  // size the partitioning vectors
    cells_per_proc.resize(nodes());
    procs_per_cell.resize(num_cells);

  // calculate the total capacity of all processors
    int total_capacity = sinit.get_capacity() * nodes();
    assert (total_capacity >= num_cells);

  // do full replication / full DD / rep-DD
    if (num_cells <= sinit.get_capacity())
    {
      // do full replication as we can fit the whole mesh on each processor
	for (int proc = 0; proc < nodes(); proc++)
	{
	    for (int cell = 1; cell <= num_cells; cell++)
	    {
		procs_per_cell[cell-1].push_back(proc);
		cells_per_proc[proc].push_back(cell);
	    }
	    assert (cells_per_proc[proc].size() == num_cells);
	}
    }
    else if (num_cells == total_capacity)
    {
      // do full DD
	int cell = 1;
	for (int proc = 0; proc < nodes(); proc++)
	{
	    int max_cell = 0;
	    while (++max_cell <= sinit.get_capacity())
	    {
		cells_per_proc[proc].push_back(cell++);
		procs_per_cell[cell-1].push_back(proc);
		assert (procs_per_cell[cell-1].size() == 1);
	    }
	    assert (procs_per_cell[proc].size() <= sinit.get_capacity());
	}
    }
    else
    {
      // do DD/replication

      // set up variables for cell replication
	int *cellrep = new int(num_cells);
	double sfrac = 0;
	int total_rep = 0;

      // loop for first estimate of cell replication
	for (int cell = 1; cell <= num_cells; cell++)
	{
	  // calculate the source fraction
	    sfrac = static_cast<double>(sinit.get_ncen(cell) +
					sinit.get_nvol(cell) +
					sinit.get_nss(cell)) / 
		(sinit.get_ncentot() + sinit.get_nsstot() + 
		 sinit.get_nvoltot());

	  // estimate the number of cell replicates where the number must be
	  // >=1 and <=number of processors
	    cellrep[cell-1] = min(nodes(), max(1, static_cast<int>
					       (sfrac * total_capacity)));
	    total_rep += cellrep[cell-1];
	}

      // calculate amount of extra capacity still available
	int xtra_reps = total_capacity - total_rep;

      // loop until xtra_capacity is gone
	while (xtra_reps != 0)
	{
	  // reset counters
	    total_rep = 0;

	    if (xtra_reps > 0)
	    {
	      // loop over cells and add one replicate at a time
		int cell = 1;
		while (cell <= num_cells && total_rep < xtra_reps)
		{
		    if (cellrep[cell-1] < nodes())
			cellrep[cell-1]++;
		    total_rep++;
		    cell++;
		}
	    }
	    else if (xtra_reps < 0)
	    {
	      // loop over cells and subtract one replicate at a time
		int cell = 1;
		while (cell <= num_cells && total_rep > xtra_reps)
		{
		    if (cellrep[cell-1] > 1)
			cellrep[cell-1]--;
		    total_rep--;
		    cell++;
		}
	    }
	    xtra_reps -= total_rep;
	}

      // now we can do our DD/replication partitioning
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    std::cout << cell << std::endl;
	    int max_rep = 0;
	    for (int proc = 0; proc < nodes(); proc++)
	    {
		std::cout << cell << std::endl;
		if (++max_rep <= cellrep[cell-1] && 
		    cells_per_proc[proc].size() < sinit.get_capacity())
		{
		    std::cout << cell << std::endl;
		    procs_per_cell[cell-1].push_back(proc);
		    std::cout << cell << std::endl;
   		    cells_per_proc[proc].push_back(cell);
		    std::cout << cell << std::endl << std::endl;
		}
 		assert (cells_per_proc[proc].size() <= sinit.get_capacity());
	    }
 	    assert (procs_per_cell[cell-1].size() > 0);
 	    assert (procs_per_cell[cell-1].size() <= nodes());
	}

      // reclaim storage
	delete [] cellrep;
    }
}

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

  // rebuilt Mesh SP
    SP<MT> return_mesh;

  // get coordinate system
    SP<Coord_sys> coord = recv_Coord();

  // get Layout
    Layout layout = recv_Layout();

  // get vertices and cell_pair
    typename MT::CCVF_a vertex = recv_vertex();
    typename MT::CCVF_i cell_pair = recv_cellpair();

  // build mesh
    return_mesh = new MT(coord, layout, vertex, cell_pair);
    std::cout << "Built MESH" << std::endl;

  // return mesh
    return return_mesh;
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

//---------------------------------------------------------------------------//
// Opacity passing interface
//---------------------------------------------------------------------------//
// send the Opacity object

template<class MT>
void Parallel_Builder<MT>::send_Opacity(const Opacity<MT> &opacity)
{
  // send out the Opacities, one component at a time

  // assure that we are on the host node
    Check (!node());

  // determine the number of cells
    int num_cells = opacity.num_cells();

  // assign the Opacity data
    double *sigma  = new double[num_cells];
    double *planck = new double[num_cells];
    double *fleck  = new double[num_cells];
    
    for (int cell = 1; cell <= num_cells; cell++)
    {
	sigma[cell-1]  = opacity.get_sigma(cell);
	planck[cell-1] = opacity.get_planck(cell);
	fleck[cell-1]  = opacity.get_fleck(cell);
    }

  // send the Opacity data
    for (int np = 1; np < nodes(); np++)
    {
	Send (num_cells, np, 20);
	Send (sigma, num_cells, np, 21);
	Send (planck, num_cells, np, 22);
	Send (fleck, num_cells, np, 23);
    }

  // delete dynamic allocation
    delete [] sigma;
    delete [] planck;
    delete [] fleck;
}

//---------------------------------------------------------------------------//
// receive the Opacity object

template<class MT>
SP<Opacity<MT> > Parallel_Builder<MT>::recv_Opacity(SP<MT> mesh)
{
  // receive and rebuild the Opacity object

  // assure we are on receive nodes
    Check (node());

  // declare return opacity object
    SP< Opacity<MT> > return_opacity;

  // check to make sure we have a valid mesh pointer
    Check (mesh);

  // receive the size of this guy
    int num_cells;
    Recv (num_cells, 0, 20);

  // check to make sure our meshes are of proper size
    Check (num_cells == mesh->num_cells());

  // make new Opacity objects
    typename MT::CCSF_double sigma(mesh);
    typename MT::CCSF_double planck(mesh);
    typename MT::CCSF_double fleck(mesh);

  // receive data from host
    double *rsigma  = new double[num_cells];
    double *rplanck = new double[num_cells];
    double *rfleck  = new double[num_cells];
    Recv (rsigma, num_cells, 0, 21);
    Recv (rplanck, num_cells, 0, 22);
    Recv (rfleck, num_cells, 0, 23);

  // assign to new Opacity objects
    for (int cell = 1; cell <= num_cells; cell++)
    {
	sigma(cell)  = rsigma[cell-1];
	planck(cell) = rplanck[cell-1];
	fleck(cell)  = rfleck[cell-1];
    }

  // reclaim dynamic memory
    delete [] rsigma;
    delete [] rplanck;
    delete [] rfleck;

  // build and return new opacity object
    return_opacity = new Opacity<MT>(sigma);
    return return_opacity;
}
    
CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Builder.cc
//---------------------------------------------------------------------------//
