//----------------------------------*-C++-*----------------------------------//
// AMR_Builder.cc
// Thomas M. Evans
// Sat Jul 25 14:25:16 1998
//---------------------------------------------------------------------------//
// @> AMR_Builder class implementation file
//---------------------------------------------------------------------------//

#include "imc/AMR_Builder.hh"
#include "imc/XYZCoord_sys.hh"
#include "imc/XYCoord_sys.hh"
#include <cmath>
#include <iostream>

IMCSPACE

// stl components 
using std::pow;
using std::endl;

//---------------------------------------------------------------------------//
// explicit constructor
//---------------------------------------------------------------------------//
// constructor (explicit)

AMR_Builder::AMR_Builder(SP<AMR_Interface> interface)
{
    Require (interface);
    
  // get data from the host
    const int *layout        = interface->get_layout();
    const double *node_coord = interface->get_node_coord();
    int num_cells            = interface->get_num_cells();

  // hardwired coordinate system
    coord_sys = "xyz";
    int dimension    = coord_sys.size();

  // assign this data to our vectors
    int laysize = num_cells * dimension * 2;
    int vsize   = num_cells * dimension * pow(static_cast<float>(2), 
					      dimension);
    laydata.resize(laysize);
    vertices.resize(vsize);
    for (int i = 0; i < laysize; i++)
	laydata[i] = layout[i];
    for (int i = 0; i < vsize; i++)
	vertices[i] = node_coord[i];
}

//---------------------------------------------------------------------------//
// build Mesh interface
//---------------------------------------------------------------------------//

SP<OS_Mesh> AMR_Builder::build_Mesh()
{
  // declare SPs
    SP<Coord_sys> coord;
    SP<Layout>    layout;
    SP<OS_Mesh>   mesh;

  // build mesh-independent objects
    coord  = build_Coord();
    layout = build_Layout(coord);
    mesh   = build_Mesh(coord, layout);

  // return mesh
    Ensure (layout->num_cells() == mesh->num_cells());
    return mesh;
}

//---------------------------------------------------------------------------//
// build mesh implementation
//---------------------------------------------------------------------------//
// build a Coord_sys

SP<Coord_sys> AMR_Builder::build_Coord()
{
  // build coordinate system
    SP<Coord_sys> coord;
    if (coord_sys == "xy" || coord_sys == "XY")
    {
	SP<XYCoord_sys> xycoord = new XYCoord_sys;
	coord = xycoord;
    }
    else if (coord_sys == "xyz" || coord_sys == "XYZ")
    {
	SP<XYZCoord_sys> xyzcoord = new XYZCoord_sys;
	coord = xyzcoord;
    }
    else
	Check (0);

  // return new coord sys
    Ensure (coord->get_Coord().size() == coord_sys.size());
    return coord;
}

//---------------------------------------------------------------------------//
// build a Layout

SP<Layout> AMR_Builder::build_Layout(SP<Coord_sys> coord)
{
  // set size of new Layout
    int dimension = coord->get_dim();
    int num_cells = laydata.size() / (dimension * 2);
    SP<Layout> layout = new Layout(num_cells);
    Layout &temp = *layout;

  // build the new layout from interface data
    int index = 0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
      // size number of faces
	layout->set_size(cell, dimension * 2);

      // assign faces for this cell
	if (dimension == 2)
	{
	    temp(cell, 1) = laydata[index];
	    temp(cell, 2) = laydata[index + dimension];
	    temp(cell, 3) = laydata[index + 1];
	    temp(cell, 4) = laydata[index + 1 + dimension];
	}
	else if (dimension == 3)
	{
	    temp(cell, 1) = laydata[index];
	    temp(cell, 2) = laydata[index + dimension];
	    temp(cell, 3) = laydata[index + 1];
	    temp(cell, 4) = laydata[index + 1 + dimension];
	    temp(cell, 5) = laydata[index + 2];
	    temp(cell, 6) = laydata[index + 2 + dimension];
	}
	else
	    Check (0);

      // update index
	index += dimension * 2;
    }
    Check (index == laydata.size());
    
  // return layout 
    return layout;
}

//---------------------------------------------------------------------------//
// build the Mesh

SP<OS_Mesh> AMR_Builder::build_Mesh(SP<Coord_sys> coord, SP<Layout> layout)
{
  // define necessary Mesh components
    int dimension = coord->get_dim();
    int num_cells = layout->num_cells();
    int num_vert  = pow(static_cast<float>(2), dimension);
    OS_Mesh::CCVF_d vertex(dimension);
    OS_Mesh::CCVF_i cell_pair(num_cells);

  // assign to vertex and cell_pair	
    int index = 0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	for (int v = 0; v < num_vert; v++)
	{
	    for (int d = 0; d < dimension; d++)
		vertex[d].push_back(vertices[index++]);
	    cell_pair[cell-1].push_back(vertex[0].size());
	}
	Check (cell_pair[cell-1].size() == num_vert);
    }
    Check (index == vertices.size());

  // build new Mesh
    SP<OS_Mesh> mesh = new OS_Mesh(coord, *layout, vertex, cell_pair);
    return mesh;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of AMR_Builder.cc
//---------------------------------------------------------------------------//
