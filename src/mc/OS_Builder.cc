//----------------------------------*-C++-*----------------------------------//
// OS_Builder.cc
// Thomas M. Evans
// Mon Feb  9 16:16:07 1998
//---------------------------------------------------------------------------//
// @> OS_Builder class implementation file
//---------------------------------------------------------------------------//

#include "OS_Builder.hh"
#include "XYCoord_sys.hh"
#include "XYZCoord_sys.hh"
#include "ds++/Assert.hh"
#include <iostream>

namespace rtt_mc 
{

using std::endl;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public Mesh build member functions
//---------------------------------------------------------------------------//

SP<OS_Mesh> OS_Builder::build_Mesh()
{
  // declare smart pointers
    SP<Coord_sys> coord;
    SP<Layout> layout;
    SP<OS_Mesh> return_mesh;
    
  // build mesh-independent objects
    coord  = build_Coord();
    layout = build_Layout(*coord);
    
  // build mesh
    int dim = coord->get_dim();
    if (dim == 2)
	return_mesh = build_2DMesh(coord, *layout);
    else if (dim == 3)
      	return_mesh = build_3DMesh(coord, *layout);
    return return_mesh;
}

//---------------------------------------------------------------------------//
// private Mesh build member functions
//---------------------------------------------------------------------------//

SP<OS_Mesh> OS_Builder::build_2DMesh(SP<Coord_sys> coord, Layout &layout)
{
  // variable declarations
    int num_xsur   = fine_edge[0].size();
    int num_ysur   = fine_edge[1].size();
    int num_xcells = num_xsur - 1;
    int num_ycells = num_ysur - 1;
    int num_cells  = num_xcells * num_ycells;
    int num_vert   = num_xsur * num_ysur;
    int dimension  = coord->get_dim();

  // check some assertions
    Check (layout.num_cells() == num_cells);

  // initialization variables for Mesh
    OS_Mesh::CCVF_d vertex(dimension);
    OS_Mesh::CCVF_i cell_pair(num_cells);

  // size vertex and cell_pair arrays, 4 vertices per cell
    for (int d = 1; d <= dimension; d++)
	vertex[d-1].resize(num_vert);
    for (int cell = 1; cell <= num_cells; cell++)
	cell_pair[cell-1].resize(4);

  // set vertex arrays
    for (int j = 1; j <= num_ysur; j++)
	for (int i = 1; i <= num_xsur; i++)
	{
	  // calculate vertex index
	    int index = 1 + (i-1) + num_xsur*(j-1);

	  // assign vertices
	    vertex[0][index-1] = fine_edge[0][i-1];
	    vertex[1][index-1] = fine_edge[1][j-1];
	}

  // set cell-pairings to vertices
    for (int j = 1; j <= num_ycells; j++)
	for (int i = 1; i <= num_xcells; i++)
	{
	  // indices for cell and lower-left vertex
	    int cell       = 1 + (i-1) + num_xcells*(j-1);
	    int ref_vertex = 1 + (i-1) + num_xsur*(j-1);

	  // pair cells to vertex indices (switch to accomodate graphics dump)
	    cell_pair[cell-1][0] = ref_vertex;
	    cell_pair[cell-1][1] = ref_vertex + 1;
	    cell_pair[cell-1][2] = ref_vertex + 1 + num_xsur;
	    cell_pair[cell-1][3] = ref_vertex + num_xsur;
	}

  // create mesh
    SP<OS_Mesh> mesh_return(new OS_Mesh(coord, layout, vertex, cell_pair));

  // return mesh to builder
    return mesh_return;
}

//---------------------------------------------------------------------------//

SP<OS_Mesh> OS_Builder::build_3DMesh(SP<Coord_sys> coord, Layout &layout)
{
  // variable declarations
    int num_xsur   = fine_edge[0].size();
    int num_ysur   = fine_edge[1].size();
    int num_zsur   = fine_edge[2].size();
    int num_xcells = num_xsur - 1;
    int num_ycells = num_ysur - 1;
    int num_zcells = num_zsur - 1;
    int num_cells  = num_xcells * num_ycells * num_zcells;
    int num_vert   = num_xsur * num_ysur * num_zsur;
    int dimension  = coord->get_dim();

  // check some assertions
    Check (layout.num_cells() == num_cells);

  // initialization variables for Mesh
    OS_Mesh::CCVF_d vertex(dimension);
    OS_Mesh::CCVF_i cell_pair(num_cells);

  // size vertex and cell_pair arrays, 8 vertices per cell
    for (int d = 1; d <= dimension; d++)
	vertex[d-1].resize(num_vert);
    for (int cell = 1; cell <= num_cells; cell++)
	cell_pair[cell-1].resize(8);

  // set vertex arrays
    for (int k = 1; k <= num_zsur; k++)
	for (int j = 1; j <= num_ysur; j++)
	    for (int i = 1; i <= num_xsur; i++)
	    {
	      // calculate vertex index
		int index = 1 + (i-1) + num_xsur*(j-1) + 
		    num_xsur*num_ysur*(k-1);

	      // assign vertices
		vertex[0][index-1] = fine_edge[0][i-1];
		vertex[1][index-1] = fine_edge[1][j-1];
		vertex[2][index-1] = fine_edge[2][k-1];
	    }

  // set cell-pairings to vertices
    for (int k = 1; k <= num_zcells; k++)
	for (int j = 1; j <= num_ycells; j++)
	    for (int i = 1; i <= num_xcells; i++)
	    {
	      // indices to cell and lower-left vertex
		int cell = 1 + (i-1) + num_xcells*(j-1) + 
		    num_xcells*num_ycells*(k-1);
		int ref_vertex = 1 + (i-1) + num_xsur*(j-1) +
		    num_xsur*num_ysur*(k-1);

	      // pair cells to vertex indices
		cell_pair[cell-1][0] = ref_vertex;
		cell_pair[cell-1][1] = ref_vertex + 1;
		cell_pair[cell-1][2] = ref_vertex + 1 + num_xsur;
		cell_pair[cell-1][3] = ref_vertex + num_xsur;
		cell_pair[cell-1][4] = cell_pair[cell-1][0] + num_xsur *
		    num_ysur;
		cell_pair[cell-1][5] = cell_pair[cell-1][1] + num_xsur *
		    num_ysur;
		cell_pair[cell-1][6] = cell_pair[cell-1][2] + num_xsur *
		    num_ysur;
		cell_pair[cell-1][7] = cell_pair[cell-1][3] + num_xsur *
		    num_ysur;
	    }

  // create mesh
    SP<OS_Mesh> mesh_return(new OS_Mesh(coord, layout, vertex, cell_pair));

  // return mesh to builder
    return mesh_return;
}

//---------------------------------------------------------------------------//
// Coord_sys build member functions
//---------------------------------------------------------------------------//

SP<Coord_sys> OS_Builder::build_Coord()
{
  // build coordinate system
    SP<Coord_sys> coord;
    if (coord_system == "xy" || coord_system == "XY")
    {
	SP<XYCoord_sys> xycoord(new XYCoord_sys);
	coord = xycoord;
    }
    else if (coord_system == "xyz" || coord_system == "XYZ")
    {
	SP<XYZCoord_sys> xyzcoord(new XYZCoord_sys);
	coord = xyzcoord;
    }

  // return base class SP to a derived Coord_sys
    return coord;
}

//---------------------------------------------------------------------------//
// Layout build member functions
//---------------------------------------------------------------------------//

SP<Layout> OS_Builder::build_Layout(const Coord_sys &coord)
{
  // set size of new Layout
    int size = 1;
    for (int d = 0; d < coord.get_dim(); d++)
	size *= fine_edge[d].size() - 1;
    SP<Layout> layout(new Layout(size));

  // set number of faces for each cell in Layout, for OS Meshes this is two
  // times the dimension of the Mesh, ie. a 2D mesh cell has 4 faces
    for (int i = 1; i <= size; i++)
	layout->set_size(i, coord.get_dim()*2);

  // assign cells and faces to Layout
    if (coord.get_dim() == 2)
	assign2D(*layout);
    else if (coord.get_dim() == 3)
	assign3D(*layout);

  // return built Layout
    return layout;
}

//---------------------------------------------------------------------------//

void OS_Builder::assign2D(Layout &layout)
{
  // 2D map of Mesh
    int num_xcells = fine_edge[0].size() - 1;
    int num_ycells = fine_edge[1].size() - 1;

  // loop over num_cells and assign cell across faces
  // 1:x(-), 2:x(+), 3:y(-), 4:y(+)
    for (int cell = 1; cell <= layout.num_cells(); cell++)
    {
	layout(cell, 1) = cell - 1;
	layout(cell, 2) = cell + 1;
	layout(cell, 3) = cell - num_xcells;
	layout(cell, 4) = cell + num_xcells;
    }

  // take care of boundary conditions
    int bcell = 0;

  // low x boundary, i = 1
    for (int j = 1; j <= num_ycells; j++)
    {
	bcell = 1 + num_xcells * (j - 1);
	if (bnd_cond[0] == "vacuum")
	    layout(bcell, 1) = 0;
	else if (bnd_cond[0] == "reflect")
	    layout(bcell, 1) = bcell;
    }

  // high x boundary, i = num_xcells
    for (int j = 1; j <= num_ycells; j++)
    {
	bcell = 1 + (num_xcells - 1) + num_xcells * (j - 1);
	if (bnd_cond[1] == "vacuum")
	    layout(bcell, 2) = 0;
	else if (bnd_cond[1] == "reflect")
	    layout(bcell, 2) = bcell;
    }

  // low y boundary, j = 1
    for (int i = 1; i <= num_xcells; i++)
    {
	bcell = 1 + (i - 1);
	if (bnd_cond[2] == "vacuum")
	    layout(bcell, 3) = 0;
	else if (bnd_cond[2] == "reflect")
	    layout(bcell, 3) = bcell;
    }

  // high y boundary, j = num_ycells
    for (int i = 1; i <= num_xcells; i++)
    {
	bcell = 1 + (i - 1) + num_xcells * (num_ycells - 1);
	if (bnd_cond[3] == "vacuum")
	    layout(bcell, 4) = 0;
	else if (bnd_cond[3] == "reflect")
	    layout(bcell, 4) = bcell;
    }
}

//---------------------------------------------------------------------------//

void OS_Builder::assign3D(Layout &layout)
{
  // 3D map of Mesh
    int num_xcells = fine_edge[0].size() - 1;
    int num_ycells = fine_edge[1].size() - 1;
    int num_zcells = fine_edge[2].size() - 1;

  // loop over num_cells and assign cell across faces
  // 1:x(-), 2:x(+), 3:y(-), 4:y(+), 5:z(-), 6:z(+)
    for (int cell = 1; cell <= layout.num_cells(); cell++)
    {
	layout(cell, 1) = cell - 1;
	layout(cell, 2) = cell + 1;
	layout(cell, 3) = cell - num_xcells;
	layout(cell, 4) = cell + num_xcells;
	layout(cell, 5) = cell - num_xcells * num_ycells;
	layout(cell, 6) = cell + num_xcells * num_ycells;
    }

  // take care of boundary conditions
    int bcell = 0;

  // low x boundary, i = 1
    for (int k = 1; k <= num_zcells; k++)
	for (int j = 1; j <= num_ycells; j++)
	{
	    bcell = 1 + num_xcells * (j - 1) + num_xcells * num_ycells * 
		(k - 1);
	    if (bnd_cond[0] == "vacuum")
		layout(bcell, 1) = 0;
	    else if (bnd_cond[0] == "reflect")
		layout(bcell, 1) = bcell;
	}

  // high x boundary, i = num_xcells
    for (int k = 1; k <= num_zcells; k++)
	for (int j = 1; j <= num_ycells; j++)
	{
	    bcell = 1 + (num_xcells - 1) + num_xcells * (j - 1) + num_xcells 
		* num_ycells * (k - 1);
	    if (bnd_cond[1] == "vacuum")
		layout(bcell, 2) = 0;
	    else if (bnd_cond[1] == "reflect")
		layout(bcell, 2) = bcell;
	}

  // low y boundary, j = 1
    for (int k = 1; k <= num_zcells; k++)
	for (int i = 1; i <= num_xcells; i++)
	{
	    bcell = 1 + (i - 1) + num_xcells * num_ycells * (k - 1);
	    if (bnd_cond[2] == "vacuum")
		layout(bcell, 3) = 0;
	    else if (bnd_cond[2] == "reflect")
		layout(bcell, 3) = bcell;
	}

  // high y boundary, j = num_ycells
    for (int k = 1; k <= num_zcells; k++)
	for (int i = 1; i <= num_xcells; i++)
	{
	    bcell = 1 + (i - 1) + num_xcells * (num_ycells - 1) +
		num_xcells * num_ycells * (k - 1);
	    if (bnd_cond[3] == "vacuum")
		layout(bcell, 4) = 0;
	    else if (bnd_cond[3] == "reflect")
		layout(bcell, 4) = bcell;
	}

  // low z boundary, k = 1
    for (int j = 1; j <= num_ycells; j++)
	for (int i = 1; i <= num_xcells; i++)
	{
	    bcell = 1 + (i - 1) + num_xcells * (j - 1);
	    if (bnd_cond[4] == "vacuum")
		layout(bcell, 5) = 0;
	    else if (bnd_cond[4] == "reflect")
		layout(bcell, 5) = bcell;
	}

  // high z boundary, k = num_zcells
    for (int j = 1; j <= num_ycells; j++)
	for (int i = 1; i <= num_xcells; i++)
	{
	    bcell = 1 + (i - 1) + num_xcells * (j - 1) + num_xcells * 
		num_ycells * (num_zcells - 1);
	    if (bnd_cond[5] == "vacuum")
		layout(bcell, 6) = 0;
	    else if (bnd_cond[5] == "reflect")
		layout(bcell, 6) = bcell;
	}
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of OS_Builder.cc
//---------------------------------------------------------------------------//
