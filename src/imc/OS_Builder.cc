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
#include <cassert>

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public Mesh build member functions
//---------------------------------------------------------------------------//
SP<OS_Mesh> OS_Builder::Build_Mesh()
{
  // declare smart pointers
    SP<Coord_sys> coord;
    SP<Layout> layout;
    SP<OS_Mesh> return_mesh;

  // get data arrays from OS_Parser needed to build OS_Mesh
    fine_edge = parser->Fine_edge();
    bnd_cond  = parser->Boundaries();
    
  // build mesh-independent objects
    coord  = Build_Coord();
    layout = Build_Layout(*coord);
    
  // build mesh
    int dim = coord->Get_dim();
    if (dim == 2)
	return_mesh = Build_2DMesh(coord, *layout);
    else if (dim == 3)
      	return_mesh = Build_3DMesh(coord, *layout);
    return return_mesh;
}

//---------------------------------------------------------------------------//
// private Mesh build member functions
//---------------------------------------------------------------------------//
SP<OS_Mesh> OS_Builder::Build_2DMesh(SP<Coord_sys> coord, Layout &layout)
{
  // variable declarations
    int num_cells  = layout.Num_cells();
    int num_xsur   = fine_edge[0].size();
    int num_ysur   = fine_edge[1].size();
    int num_xcells = num_xsur - 1;
    int num_ycells = num_ysur - 1;
    int cell       = 0;
    int dimension  = coord->Get_dim();

  // initialization variables for Mesh
    OS_Mesh::CCVF_a pos(2);
    OS_Mesh::CCVF_a dim(2);

  // size position and dimension arrays
    for (int d = 1; d <= dimension; d++)
    {
	pos[d-1].resize(num_cells);
	dim[d-1].resize(num_cells);
    }

  // set position and dimension arrays
    for (int i = 1; i <= num_xcells; i++)
	for (int j = 1; j <= num_ycells; j++)
	{
	    cell           = 1 + (i-1) + num_xcells * (j-1);
	    pos[0][cell-1] = (fine_edge[0][i-1] + fine_edge[0][i]) / 2.0;
	    pos[1][cell-1] = (fine_edge[1][j-1] + fine_edge[1][j]) / 2.0;
	    dim[0][cell-1] = fine_edge[0][i] - fine_edge[0][i-1];
	    dim[1][cell-1] = fine_edge[1][j] - fine_edge[1][j-1];
	}

  // return mesh to builder
    SP<OS_Mesh> mesh_return = new OS_Mesh(coord, layout, pos, dim, fine_edge);
    return mesh_return;
}

SP<OS_Mesh> OS_Builder::Build_3DMesh(SP<Coord_sys> coord, Layout &layout)
{
  // variable declarations
    int num_cells  = layout.Num_cells();
    int num_xsur   = fine_edge[0].size();
    int num_ysur   = fine_edge[1].size();
    int num_zsur   = fine_edge[2].size();
    int num_xcells = num_xsur - 1;
    int num_ycells = num_ysur - 1;
    int num_zcells = num_zsur - 1;
    int cell       = 0;
    int dimension  = coord->Get_dim();

  // initialization variables for Mesh
    OS_Mesh::CCVF_a pos(3);
    OS_Mesh::CCVF_a dim(3);

  // size position and dimension arrays
    for (int d = 1; d <= dimension; d++)
    {
	pos[d-1].resize(num_cells);
	dim[d-1].resize(num_cells);
    }

  // set position and dimension arrays
    for (int i = 1; i <= num_xcells; i++)
	for (int j = 1; j <= num_ycells; j++)
	    for (int k = 1; k <= num_zcells; k++)
	    {
		cell = 1 + (i-1) + num_xcells * (j-1) + num_xcells * 
		    num_ycells * (k-1);
		pos[0][cell-1] = (fine_edge[0][i-1] + fine_edge[0][i]) / 2.0;
		pos[1][cell-1] = (fine_edge[1][j-1] + fine_edge[1][j]) / 2.0;
		pos[2][cell-1] = (fine_edge[2][k-1] + fine_edge[2][k]) / 2.0;
		dim[0][cell-1] = fine_edge[0][i] - fine_edge[0][i-1];
		dim[1][cell-1] = fine_edge[1][j] - fine_edge[1][j-1];
		dim[2][cell-1] = fine_edge[2][k] - fine_edge[2][k-1];
	    }

  // return mesh to builder
    SP<OS_Mesh> mesh_return = new OS_Mesh(coord, layout, pos, dim, fine_edge);
    return mesh_return;
}

//---------------------------------------------------------------------------//
// Coord_sys build member functions
//---------------------------------------------------------------------------//
SP<Coord_sys> OS_Builder::Build_Coord()
{
  // get coordinate system from OS_Parser
    string coord_system = parser->Coordinates();

  // build coordinate system
    SP<Coord_sys> coord;
    if (coord_system == "xy" || coord_system == "XY")
    {
	SP<XYCoord_sys> xycoord = new XYCoord_sys;
	coord = xycoord;
    }
    else if (coord_system == "xyz" || coord_system == "XYZ")
    {
	SP<XYZCoord_sys> xyzcoord = new XYZCoord_sys;
	coord = xyzcoord;
    }

  // return base class SP to a derived Coord_sys
    return coord;
}

//---------------------------------------------------------------------------//
// Layout build member functions
//---------------------------------------------------------------------------//
SP<Layout> OS_Builder::Build_Layout(const Coord_sys &coord)
{
  // set size of new Layout
    int size = 1;
    for (int d = 0; d < coord.Get_dim(); d++)
	size *= fine_edge[d].size() - 1;
    SP<Layout> layout = new Layout(size);

  // set number of faces for each cell in Layout, for OS Meshes this is two
  // times the dimension of the Mesh, ie. a 2D mesh cell has 4 faces
    for (int i = 1; i < size; i++)
	layout->Set_size(i, coord.Get_dim()*2);

  // assign cells and faces to Layout
    if (coord.Get_dim() == 2)
	Assign2D(*layout);
    else if (coord.Get_dim() == 3)
	Assign3D(*layout);

  // return built Layout
    return layout;
}

void OS_Builder::Assign2D(Layout &layout)
{
  // 2D map of Mesh
    int num_xcells = fine_edge[0].size() - 1;
    int num_ycells = fine_edge[1].size() - 1;

  // loop over num_cells and assign cell across faces
  // 1:x(-), 2:x(+), 3:y(-), 4:y(+)
    for (int cell = 1; cell <= layout.Num_cells(); cell++)
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

void OS_Builder::Assign3D(Layout &layout)
{
  // 3D map of Mesh
    int num_xcells = fine_edge[0].size() - 1;
    int num_ycells = fine_edge[1].size() - 1;
    int num_zcells = fine_edge[2].size() - 1;

  // loop over num_cells and assign cell across faces
  // 1:x(-), 2:x(+), 3:y(-), 4:y(+), 5:z(-), 6:z(+)
    for (int cell = 1; cell <= layout.Num_cells(); cell++)
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

CSPACE

//---------------------------------------------------------------------------//
//                              end of OS_Builder.cc
//---------------------------------------------------------------------------//
