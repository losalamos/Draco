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
    
  // call parser
    parser();
    
  // build mesh-independent objects
    coord  = buildCoord();
    layout = buildLayout(*coord);
    
  // build mesh
    int dim = coord->getDim();
    if (dim == 2)
	return_mesh = build2DMesh(coord, *layout);
    else if (dim == 3)
      //	return_mesh = build3DMesh(coord, *layout);
	return_mesh = 0;
    return return_mesh;
}

//---------------------------------------------------------------------------//
// private Mesh build member functions
//---------------------------------------------------------------------------//
SP<OS_Mesh> OS_Builder::Build_2DMesh(SP<Coord_sys> coord, 
				     const Layout &layout)
{
  // variable declarations
    int num_cells  = layout.getNum_cell();
    int num_xsur   = fine_edge[0].size();
    int num_ysur   = fine_edge[1].size();
    int num_xcells = num_xsur - 1;
    int num_ycells = num_ysur - 1;
    int cell       = 0;
    int dimension  = coord->getDim();

  // initialization variables for Mesh
    OS_Mesh::CCVF_a pos(2);
    OS_Mesh::CCVF_a dim(2);
    OS_Mesh::CCF_i  index(num_cells);

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
	    index[cell-1]  = cell;
	}

  // return mesh to builder
    SP<OS_Mesh> mesh_return = new OS_Mesh(coord, layout, pos, dim, fine_edge,
					  index);
    return mesh_return;
}

//---------------------------------------------------------------------------//
// parser member functions
//---------------------------------------------------------------------------//
void OS_Builder::Parser()
{
  // open input file, ifstream object requires C-style string
    const char *file = input_file.c_str();
    ifstream input(file);
    
  // determine coord_sys
    string keyword;
    while (keyword != "end-title")
    {
      // test that we have not reached end-of-file
	assert (!input.eof());

	input >> keyword;
	if (keyword == "coord:")
	    input >> coord_system;
    }

  // check to make sure an appropriate coord_sys is called
    assert (coord_system == "xy" || coord_system == "XY" || 
	    coord_system == "rz" || coord_system == "RZ" ||
	    coord_system == "XYZ" || coord_system == "xyz");

  // call appropriate sub-parser
    if (coord_system == "xy" || coord_system == "XY" || coord_system == "rz"
	|| coord_system == "RZ")
    {
	fine_cells.resize(2);
	coarse_edge.resize(2);
	fine_edge.resize(2);
	bnd_cond.resize(4);
	parser2D(input);
    }
    else if (coord_system == "xyz" || coord_system == "XYZ")
    {
	fine_cells.resize(3);
	coarse_edge.resize(3);
	fine_edge.resize(3);
	bnd_cond.resize(6);
	parser3D(input);
    }

  // calculate fine_edge array

  // determine size of fine_edge[i] arrays
    for (int d = 0; d < fine_edge.size(); d++)
    {
	int nfine = 0;
	for (int n = 0; n < fine_cells[d].size(); n++)
	    nfine += fine_cells[d][n];
	fine_edge[d].resize(nfine+1);
    }

  // assign edge data to arrays
    for (int d = 0; d < fine_edge.size(); d++)
    {
	int ifine    = 0;
	double delta = 0.0;
	for (int i = 0; i < coarse_edge[d].size() - 1; i++)
	{
	    delta = (coarse_edge[d][i+1] - coarse_edge[d][i]) /
		fine_cells[d][i];
	    fine_edge[d][ifine] = coarse_edge[d][i];
	    for (int j = 1; j <= fine_cells[d][i]; j++)
	    {
		ifine++;
		fine_edge[d][ifine] = coarse_edge[d][i] + j * delta;
	    }
	}
      	fine_edge[d][ifine] = coarse_edge[d].back();
    }
}

void OS_Builder::Parser2D(ifstream &in)
{
  // 2D parser

    string keyword;
    int data;
    
  // initialization block input
    while (keyword != "end-init")
    {
      // test that we have not reached end-of-file
	assert (!in.eof());

      // do input
	in >> keyword;
	if (keyword == "num_xcoarse:" || keyword == "num_rcoarse:")
	{
	    in >> data;
	    fine_cells[0].resize(data);
	    coarse_edge[0].resize(data+1);  
	}
	if (keyword == "num_ycoarse:" || keyword == "num_zcoarse:")
	{
	    in >> data;
	    fine_cells[1].resize(data);
	    coarse_edge[1].resize(data+1);
	}
	if (keyword == "lox_bnd:" || keyword == "lor_bnd:")
	    in >> bnd_cond[0];
	if (keyword == "hix_bnd:" || keyword == "hir_bnd:")
	    in >> bnd_cond[1];
	if (keyword == "loy_bnd:" || keyword == "loz_bnd:")
	    in >> bnd_cond[2];
	if (keyword == "hiy_bnd:" || keyword == "hiz_bnd:")
	    in >> bnd_cond[3];
    }

  // mesh block input
    while (keyword != "end-mesh")
    {
      // test that we have not reached end-of-file
	assert (!in.eof());

      // do input
	in >> keyword;
	if (keyword == "xcoarse:" || keyword == "rcoarse:")
	    for (int i = 0; i < coarse_edge[0].size(); i++)
		in >> coarse_edge[0][i];
	if (keyword == "num_xfine:" || keyword == "num_rfine:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
		in >> fine_cells[0][i];
	if (keyword == "ycoarse:" || keyword == "zcoarse:")
	    for (int i = 0; i < coarse_edge[1].size(); i++)
		in >> coarse_edge[1][i];
	if (keyword == "num_yfine:" || keyword == "num_zfine:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
		in >> fine_cells[1][i];
    }
}    

void OS_Builder::Parser3D(ifstream &in)
{
  // 3D parser

    string keyword;
    int data;
    
  // initialization block input
    while (keyword != "end-init")
    {
      // test that we have not reached end-of-file
	assert (!in.eof());

      // do input
	in >> keyword;
	if (keyword == "num_xcoarse:")
	{
	    in >> data;
	    fine_cells[0].resize(data);
	    coarse_edge[0].resize(data+1);  
	}
	if (keyword == "num_ycoarse:")
	{
	    in >> data;
	    fine_cells[1].resize(data);
	    coarse_edge[1].resize(data+1);
	}
	if (keyword == "num_zcoarse:")
	{
	    in >> data;
	    fine_cells[2].resize(data);
	    coarse_edge[2].resize(data+1);
	}
	if (keyword == "lox_bnd:")
	    in >> bnd_cond[0];
	if (keyword == "hix_bnd:")
	    in >> bnd_cond[1];
	if (keyword == "loy_bnd:")
	    in >> bnd_cond[2];
	if (keyword == "hiy_bnd:")
	    in >> bnd_cond[3];
	if (keyword == "loz_bnd:")
	    in >> bnd_cond[4];
	if (keyword == "hiz_bnd:")
	    in >> bnd_cond[5];
    }

  // mesh block input
    while (keyword != "end-mesh")
    {
      // test that we have not reached end-of-file
	assert (!in.eof());
     
      // do input
	in >> keyword;
	if (keyword == "xcoarse:")
	    for (int i = 0; i < coarse_edge[0].size(); i++)
		in >> coarse_edge[0][i];
	if (keyword == "num_xfine:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
		in >> fine_cells[0][i];
	if (keyword == "ycoarse:")
	    for (int i = 0; i < coarse_edge[1].size(); i++)
		in >> coarse_edge[1][i];
	if (keyword == "num_yfine:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
		in >> fine_cells[1][i];
	if (keyword == "zcoarse:")
	    for (int i = 0; i < coarse_edge[2].size(); i++)
		in >> coarse_edge[2][i];
	if (keyword == "num_zfine:")
	    for (int i = 0; i < fine_cells[2].size(); i++)
		in >> fine_cells[2][i];
    }
}  

//---------------------------------------------------------------------------//
// Coord_sys build member functions
//---------------------------------------------------------------------------//
SP<Coord_sys> OS_Builder::Build_Coord()
{
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
    for (int d = 0; d < coord.getDim(); d++)
	size *= fine_edge[d].size() - 1;
    SP<Layout> layout = new Layout(size);

  // set number of faces for each cell in Layout, for OS Meshes this is two
  // times the dimension of the Mesh, ie. a 2D mesh cell has 4 faces
    for (int i = 1; i < size; i++)
	layout->setSize(i, coord.getDim()*2);

  // assign cells and faces to Layout
    if (coord.getDim() == 2)
	Assign2D(*layout);
    else if (coord.getDim() == 3)
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
    for (int cell = 1; cell <= layout.getNum_cell(); cell++)
    {
	layout(cell, 1) = cell - 1;
	layout(cell, 2) = cell + 1;
	layout(cell, 3) = cell - num_xcells;
	layout(cell, 4) = cell + num_xcells;
    }

  // take care of boundary conditions

  // low x boundary
    {
	int bcell = 1;
	for (int i = 1; i <= num_ycells; i++)
	{
	    if (bnd_cond[0] == "vacuum")
		layout(bcell, 1) = 0;
	    else if (bnd_cond[0] == "reflect")
		layout(bcell, 1) = bcell;
	    bcell += num_xcells;
	}
    }
  // high x boundary
    {
	int bcell = num_xcells;
	for (int i = 1; i <= num_ycells; i++)
	{
	    if (bnd_cond[1] == "vacuum")
		layout(bcell, 2) = 0;
	    else if (bnd_cond[1] == "reflect")
		layout(bcell, 2) = bcell;
	    bcell += num_xcells;
	}
    }
  // low y boundary
    {
	int bcell = 1;
	for (int i = 1; i <= num_xcells; i++)
	{
	    if (bnd_cond[2] == "vacuum")
		layout(bcell, 3) = 0;
	    else if (bnd_cond[2] == "reflect")
		layout(bcell, 3) = bcell;
	    bcell += 1;
	}
    }
  // high y boundary
    {
	int bcell = 1 + num_xcells * (num_ycells - 1);
	for (int i = 1; i <= num_xcells; i++)
	{
	    if (bnd_cond[3] == "vacuum")
		layout(bcell, 4) = 0;
	    else if (bnd_cond[3] == "reflect")
		layout(bcell, 4) = bcell;
	    bcell += 1;
	}
    }
}

void OS_Builder::Assign3D(Layout &layout)
{
  // need to add 3D assign function
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of OS_Builder.cc
//---------------------------------------------------------------------------//
