//----------------------------------*-C++-*----------------------------------//
// OS_Builder.cc
// Thomas M. Evans
// Mon Feb  9 16:16:07 1998
//---------------------------------------------------------------------------//
// @> OS_Builder class implementation file
//---------------------------------------------------------------------------//

#include "OS_Builder.hh"

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// build mesh public member function
//---------------------------------------------------------------------------//
SP<OS_Mesh> OS_Builder::buildMesh()
{
  // declare smart pointers
    SP<Coord_sys> coord;
    SP<Layout> layout;
    SP<OS_Mesh> return_mesh;

  // call parser
    parser();

  // build mesh-independent objects
    coord       = buildCoord();
    layout      = buildLayout(*coord);

  // build mesh
    int dim = coord->getDim();
    if (dim == 2)
	return_mesh = build2DMesh(coord, *layout);
    else if (dim == 3)
	return_mesh = build3DMesh(coord, *layout);
    return return_mesh;
}

//---------------------------------------------------------------------------//
// parser member functions
//---------------------------------------------------------------------------//
void OS_Builder::parser()
{
  // open input file
    const char *file = input_file.c_str();
    ifstream input(file);
    
  // determine coord_sys
    string keyword;
    input >> keyword >> coord_system;

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
}

void OS_Builder::parser2D(ifstream &in)
{
  // 2D parser

    string keyword;
    int data;
    
  // initialization block input
    while (keyword != "end-init")
    {
	in >> keyword;
	if (keyword == "xcoarse:" || keyword == "rcoarse:")
	{
	    in >> data;
	    fine_cells[0].resize(data);
	    coarse_edge[0].resize(data+1);  
	}
	if (keyword == "ycoarse:" || keyword == "zcoarse:")
	{
	    in >> data;
	    fine_cells[1].resize(data);
	    coarse_edge[1].resize(data+1);
	}
	if (keyword == "xfine:" || keyword == "rfine:")
	{
	    in >> data;
	    fine_edge[0].resize(data+1);
	}
	if (keyword == "yfine:" || keyword == "zfine:")
	{
	    in >> data;
	    fine_edge[1].resize(data+1);
	}
    }

  // mesh block input
    while (keyword != "end-mesh")
    {
	in >> keyword;
	if (keyword == "xdim:" || keyword == "rdim:")
	    for (int i = 0; i < coarse_edge[0].size(); i++)
		in >> coarse_edge[0][i];
	if (keyword == "xints:" || keyword == "rints:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
		in >> fine_cells[0][i];
	if (keyword == "ydim:" || keyword == "zdim:")
	    for (int i = 0; i < coarse_edge[1].size(); i++)
		in >> coarse_edge[1][i];
	if (keyword == "yints:" || keyword == "zints:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
		in >> fine_cells[1][i];
    }
}    

void OS_Builder::parser3D(ifstream &in)
{
  // 3D parser

    string keyword;
    int data;
    
  // initialization block input
    while (keyword != "end-init")
    {
	in >> keyword;
	if (keyword == "xcoarse:")
	{
	    in >> data;
	    fine_cells[0].resize(data);
	    coarse_edge[0].resize(data+1);  
	}
	if (keyword == "ycoarse:")
	{
	    in >> data;
	    fine_cells[1].resize(data);
	    coarse_edge[1].resize(data+1);
	}
	if (keyword == "zcoarse:")
	{
	    in >> data;
	    fine_cells[2].resize(data);
	    coarse_edge[2].resize(data+1);
	}
	if (keyword == "xfine:")
	{
	    in >> data;
	    fine_edge[0].resize(data+1);
	}
	if (keyword == "yfine:")
	{
	    in >> data;
	    fine_edge[1].resize(data+1);
	}
	if (keyword == "zfine:")
	{
	    in >> data;
	    fine_edge[2].resize(data+1);
	}
    }

  // mesh block input
    while (keyword != "end-mesh")
    {
	in >> keyword;
	if (keyword == "xdim:")
	    for (int i = 0; i < coarse_edge[0].size(); i++)
		in >> coarse_edge[0][i];
	if (keyword == "xints:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
		in >> fine_cells[0][i];
	if (keyword == "ydim:")
	    for (int i = 0; i < coarse_edge[1].size(); i++)
		in >> coarse_edge[1][i];
	if (keyword == "yints:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
		in >> fine_cells[1][i];
	if (keyword == "zdim:")
	    for (int i = 0; i < coarse_edge[2].size(); i++)
		in >> coarse_edge[2][i];
	if (keyword == "zints:")
	    for (int i = 0; i < fine_cells[2].size(); i++)
		in >> fine_cells[2][i];
    }
}  

CSPACE

//---------------------------------------------------------------------------//
//                              end of OS_Builder.cc
//---------------------------------------------------------------------------//
