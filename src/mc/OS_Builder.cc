//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/OS_Builder.cc
 * \author Thomas M. Evans
 * \date   Mon Feb  9 16:16:07 1998
 * \brief  OS_Builder implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "OS_Builder.hh"
#include "XYCoord_sys.hh"
#include "XYZCoord_sys.hh"
#include "ds++/Assert.hh"

#include <algorithm>

namespace rtt_mc 
{

using rtt_dsxx::SP;

using std::string;
using std::vector;
using std::endl;
using std::fill;
using std::fabs;

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Build an OS_Mesh from data specified in the OS_Mesh format file.
 *
 * The builder takes data from the OS_Mesh input format that was parsed in
 * the constructor and builds an OS_Mesh.  The builder checks for the
 * existence of the built mesh; thus, a builder can only build one mesh.  To
 * get extra copies of the mesh call the get_Mesh() accessor function.
 */
OS_Builder::SP_Mesh OS_Builder::build_Mesh()
{
    Require (!mesh);
    
    // declare smart pointers
    SP_Coord_sys coord;
    SP_Layout    layout;
    
    // build mesh-independent objects
    coord  = build_Coord();
    layout = build_Layout(*coord);
    
    // build mesh
    int dim = coord->get_dim();
    if (dim == 2)
	mesh = build_2DMesh(coord, *layout);
    else if (dim == 3)
	mesh = build_3DMesh(coord, *layout);
    
    // calculate defined surface cells
    calc_defined_surcells();
    
    // do some checks and return completed mesh
    Ensure (mesh->full_Mesh());
    Ensure (mesh->get_Coord().get_dim() == fine_edge.size());
    Ensure (mesh->num_cells() == zone.size());

    return mesh;
}

//---------------------------------------------------------------------------//
// Return the cell regions for graphics dumping

OS_Builder::sf_int OS_Builder::get_regions() const
{
    vector<int> return_regions(zone.size());
    vector<int>::const_iterator cell_itr;

    // if no user defined regions given, simply return region 1 in each cell
    if (regions.size() == 1)
    {
	// make all regions 1
	std::fill(return_regions.begin(), return_regions.end(), 1);
    }
    else if (regions.size() > 1)
    {
	// preliminary declarations
	int zone;
	for (int i = 0; i < regions.size(); i++)
	    for (int j = 0; j < regions[i].size(); j++)
	    {
		// get a zone in this region
		zone = regions[i][j];

		// find the cells in this zone and add them to the region
		for (int k = 0; k < cell_zone[zone-1].size(); k++)
		    return_regions[cell_zone[zone-1][k]-1] = i+1; 
	    }
    }
    else
	Insist (0, "Something wrong with the regions setting!");

    // make sure there are no undefined regions
    cell_itr = std::find(return_regions.begin(), return_regions.end(), 0);
    Insist (cell_itr == return_regions.end(), "Incomplete regions defined!");
    
    // return the regions
    return return_regions;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the list of defined surface cells.
 */
OS_Builder::vf_int OS_Builder::get_defined_surcells() const
{
    Require (mesh);
    return defined_surcells;
}

//---------------------------------------------------------------------------//
// PRIVATE MESH FILE PARSING IMPLEMENTATION
//---------------------------------------------------------------------------//
// Basic Mesh parser.

void OS_Builder::parser()
{
    // first open the input file
    const char *file = mesh_file.c_str();
    std_ifstream input(file);

    // make sure the file exists
    Insist (input, "You must supply a mesh input file!");

    // determine coord_sys
    string keyword;
    while (keyword != "end-title")
    {
	// test that we have not reached end-of-file
	Insist (!input.eof(), "You did not supply an end-title!");

	// do title block input
	input >> keyword;
	if (keyword == "coord:")
	    input >> coord_system;
    }

    // check to make sure an appropriate coord_sys is called
    Insist (coord_system == "xy"  || coord_system == "XY" || 
	    coord_system == "rz"  || coord_system == "RZ" ||
	    coord_system == "XYZ" || coord_system == "xyz", 
	    "Invalid coordinate system choosen!");

    // call appropriate sub-parser
    if (coord_system == "xy" || coord_system == "XY" || coord_system == "rz"
	|| coord_system == "RZ")
    {
	// size the mesh data
	fine_cells.resize(2);
	fine_ratio.resize(2);
	accum_cells.resize(2);
	coarse_edge.resize(2);
	fine_edge.resize(2);
	bnd_cond.resize(4);

	// the boundary cond is always reflecting at r=0 for an RZ Mesh
	if (coord_system == "rz" || coord_system == "RZ")
	    bnd_cond[0] = "reflect";
	
	// parse mesh particulars
	parser2D(input);
    }
    else if (coord_system == "xyz" || coord_system == "XYZ")
    {
	// size the mesh data
	fine_cells.resize(3);
	fine_ratio.resize(3);
	accum_cells.resize(3);
	coarse_edge.resize(3);
	fine_edge.resize(3);
	bnd_cond.resize(6);

	// parse mesh particulars
	parser3D(input);
    }

    // parse the source block that contains surface source information
    source_parser(input);

    // size regions data if not already sized
    if (regions.empty())
	regions.resize(1);

    // >>> calculate fine_edge array

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
        // initialize fine cell counter 
	int ifine = 0;

        // initialize local-use widths
	double coarse_cell_width     = 0.0;
	double first_fine_cell_width = 0.0;

        // loop over coarse cells -- i-th coarse cell bracketed by the i-th
        // and i+1st coarse cell edge 
	for (int i = 0; i < coarse_edge[d].size() - 1; i++)
	{
            // set the first fine edge -- equiv to first coarse edge
	    fine_edge[d][ifine] = coarse_edge[d][i];

	    // calculate the width of this coarse cell
	    coarse_cell_width = coarse_edge[d][i+1] - coarse_edge[d][i];

	    // set the width of the first fine cell in this coarse cell 
	    if (fine_ratio[d][i] == 1.0)
		first_fine_cell_width = coarse_cell_width / fine_cells[d][i]; 
	    else 
		first_fine_cell_width = (1.0-fine_ratio[d][i]) *
		    coarse_cell_width /
		    (1.0-pow(fine_ratio[d][i],fine_cells[d][i])); 

	    // initialize the offset from the low edge of the coarse cell
	    double coarse_offset = 0.0;

	    // set each fine edge in this coarse cell
	    for (int j = 1; j <= fine_cells[d][i]; j++)
	    {
	        // increment fine edge counter (for entire direction d)
		ifine++;

	        // set offset from low edge of this coarse cell; set fine
	        // edge 
		if (fine_ratio[d][i] == 1.0)
		    coarse_offset = j * first_fine_cell_width;
		else
		    coarse_offset += pow(fine_ratio[d][i],j-1) *
			first_fine_cell_width; 
		fine_edge[d][ifine] = coarse_edge[d][i] + coarse_offset;
	    }

	    // should really set the last fine edge in this coarse cell 
	    // to the next coarse edge, but it could affect exact 
	    // regression tracking/testing.
	    // fine_edge[d][ifine] = coarse_edge[d][i+1];

	    // check last fine cell in this coarse cell; it should be equal
	    // to the "first fine cell width" if we use the inverse of the
	    // ratio. 
	    double last_width = fine_edge[d][ifine] - fine_edge[d][ifine-1]; 
	    double inv_ratio = 1.0 / fine_ratio[d][i];
	    double check_last_w;
	    if (inv_ratio == 1.0)
		check_last_w = coarse_cell_width / fine_cells[d][i];
	    else
		check_last_w = (1.0-inv_ratio) * coarse_cell_width /
		    (1.0-pow(inv_ratio,fine_cells[d][i]));
	    Check (fabs(last_width - check_last_w) < 1.0e-5 * last_width);
	}

        // reset the last fine edge -- equiv to the last coarse edge.  
      	fine_edge[d][ifine] = coarse_edge[d].back();
    }
}

//---------------------------------------------------------------------------//
// Parse a 2D mesh.

void OS_Builder::parser2D(std_ifstream &in)
{
    // 2D parser

    string keyword;
    int data;
    
    // initialization block input
    while (keyword != "end-init")
    {
	// test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-init!");

	// do input
	in >> keyword;
	if (keyword == "num_xcoarse:" || keyword == "num_rcoarse:")
	{
	    in >> data;
	    fine_cells[0].resize(data);
	    fine_ratio[0].resize(data);
	    fill(fine_ratio[0].begin(), fine_ratio[0].end(), 1.0);
	    coarse_edge[0].resize(data+1); 
	}
	if (keyword == "num_ycoarse:" || keyword == "num_zcoarse:")
	{
	    in >> data;
	    fine_cells[1].resize(data);
	    fine_ratio[1].resize(data);
	    fill(fine_ratio[1].begin(), fine_ratio[1].end(), 1.0);
	    coarse_edge[1].resize(data+1);
	}
	if (keyword == "lox_bnd:" || keyword == "lor_bnd:")
	{
	    in >> bnd_cond[0];
	    Insist (((coord_system == "rz") || (coord_system == "RZ")) ?
		    bnd_cond[0] == "reflect" : true, "Invalid b.c. for RZ!"); 
	}
	if (keyword == "hix_bnd:" || keyword == "hir_bnd:")
	    in >> bnd_cond[1];
	if (keyword == "loy_bnd:" || keyword == "loz_bnd:")
	    in >> bnd_cond[2];
	if (keyword == "hiy_bnd:" || keyword == "hiz_bnd:")
	    in >> bnd_cond[3];
	if (keyword == "num_regions:")
	{
	    in >> data;
	    regions.resize(data);
	}
	if (keyword == "num_zones_per_region:")
	{
	    Insist (regions.size() > 0, "Region size not set!");
	    for (int i = 0; i < regions.size(); i++)
	    {
		in >> data;
		regions[i].resize(data);
	    }
	}
	if (keyword == "regions:")
	{
	    Insist (regions.size() > 0, "Region size not set!");
	    in >> data;
	    for (int i = 0; i < regions[data-1].size(); i++)
		in >> regions[data-1][i];
	}
    }

    // mesh block input
    while (keyword != "end-mesh")
    {
	// test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-mesh!");

	// do input
	in >> keyword;
	if (keyword == "xcoarse:" || keyword == "rcoarse:")
	{
	    for (int i = 0; i < coarse_edge[0].size(); i++)
		in >> coarse_edge[0][i];
	    Insist (((coord_system == "rz") || (coord_system == "RZ")) ?
		    coarse_edge[0][0] == 0.0 : true, 
		    "The first radial value must be zero!"); 
	}
	if (keyword == "num_xfine:" || keyword == "num_rfine:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
		in >> fine_cells[0][i];
	if (keyword == "xfine_ratio:" || keyword == "rfine_ratio:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
	    {
		in >> fine_ratio[0][i];
		Check (fine_ratio[0][i] > 0.0);
	    }
	if (keyword == "ycoarse:" || keyword == "zcoarse:")
	    for (int i = 0; i < coarse_edge[1].size(); i++)
		in >> coarse_edge[1][i];
	if (keyword == "num_yfine:" || keyword == "num_zfine:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
		in >> fine_cells[1][i];
	if (keyword == "yfine_ratio:" || keyword == "zfine_ratio:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
	    {
		in >> fine_ratio[1][i];
		Check (fine_ratio[1][i] > 0.0);
	    }
	
	// unfolding angle for the RZWedge_Mesh
	if (keyword == "wedge_angle_degrees:")
	{
	    Insist ((coord_system == "rz") || (coord_system == "RZ"),
		    "Wedge angle specified, but coord_sys not equal to rz!"); 
	    in >> theta_degrees;
	    Insist (theta_degrees  >  0.0, "RZWedge angle <= 0!");
	    Insist (theta_degrees <= 90.0, "RZWedge angle > 90!");
	}    
    }
} 

//---------------------------------------------------------------------------//
// Parse a 3D mesh.

void OS_Builder::parser3D(std_ifstream &in)
{
    // 3D parser

    string keyword;
    int data;
    
    // initialization block input
    while (keyword != "end-init")
    {
	// test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-init!");

	// do input
	in >> keyword;
	if (keyword == "num_xcoarse:")
	{
	    in >> data;
	    fine_cells[0].resize(data);
	    fine_ratio[0].resize(data);
	    fill(fine_ratio[0].begin(), fine_ratio[0].end(), 1.0);
	    coarse_edge[0].resize(data+1);  
	}
	if (keyword == "num_ycoarse:")
	{
	    in >> data;
	    fine_cells[1].resize(data);
	    fine_ratio[1].resize(data);
	    fill(fine_ratio[1].begin(), fine_ratio[1].end(), 1.0);
	    coarse_edge[1].resize(data+1);
	}
	if (keyword == "num_zcoarse:")
	{
	    in >> data;
	    fine_cells[2].resize(data);
	    fine_ratio[2].resize(data);
	    fill(fine_ratio[2].begin(), fine_ratio[2].end(), 1.0);
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
	if (keyword == "num_regions:")
	{
	    in >> data;
	    regions.resize(data);
	}
	if (keyword == "num_zones_per_region:")
	{
	    Insist (regions.size() > 0, "Region size not set!");
	    for (int i = 0; i < regions.size(); i++)
	    {
		in >> data;
		regions[i].resize(data);
	    }
	}
	if (keyword == "regions:")
	{
	    Insist (regions.size() > 0, "Region size not set!");
	    in >> data;
	    for (int i = 0; i < regions[data-1].size(); i++)
		in >> regions[data-1][i];
	}
    }

    // mesh block input
    while (keyword != "end-mesh")
    {
	// test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-mesh!");
     
	// do input
	in >> keyword;
	if (keyword == "xcoarse:")
	    for (int i = 0; i < coarse_edge[0].size(); i++)
		in >> coarse_edge[0][i];
	if (keyword == "num_xfine:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
		in >> fine_cells[0][i];
	if (keyword == "xfine_ratio:")
	    for (int i = 0; i < fine_cells[0].size(); i++)
	    {
		in >> fine_ratio[0][i];
		Check (fine_ratio[0][i] > 0.0);
	    }
	if (keyword == "ycoarse:")
	    for (int i = 0; i < coarse_edge[1].size(); i++)
		in >> coarse_edge[1][i];
	if (keyword == "num_yfine:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
		in >> fine_cells[1][i];
	if (keyword == "yfine_ratio:")
	    for (int i = 0; i < fine_cells[1].size(); i++)
	    {
		in >> fine_ratio[1][i];
		Check (fine_ratio[1][i] > 0.0);
	    }
	if (keyword == "zcoarse:")
	    for (int i = 0; i < coarse_edge[2].size(); i++)
		in >> coarse_edge[2][i];
	if (keyword == "num_zfine:")
	    for (int i = 0; i < fine_cells[2].size(); i++)
		in >> fine_cells[2][i];
	if (keyword == "zfine_ratio:")
	    for (int i = 0; i < fine_cells[2].size(); i++)
	    {
		in >> fine_ratio[2][i];
		Check (fine_ratio[2][i] > 0.0);
	    }
    }
} 

//---------------------------------------------------------------------------//
// Parse source descriptions that are part of the mesh format.

void OS_Builder::source_parser(std_ifstream &in)
{
    // input keywords
    string keyword;
    int    data;
  
    // read source info, at present this information is exclusively surface
    // source position information
    while (keyword != "end-source")
    {
	// test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-source!");

	// do input
	in >> keyword;

	// check keyword and read appropriate data

	// get the number of surface sources
	if (keyword == "num_ss:")
	{
	    in >> data;
	    
	    // size surface source data
	    ss_pos.resize(data);
	    num_defined_surcells.resize(data);
	    defined_surcells.resize(data);
	    
	    // initialize to zero
	    fill(num_defined_surcells.begin(), 
		 num_defined_surcells.end(), 0); 
	}

	// read in the surface source positions
	if (keyword == "sur_source:")
	{
	    // make sure that surface source has been sized
	    Insist (ss_pos.size() != 0, "Number of surface sources = 0!"); 

	    // input the surface source position
	    for (int i = 0; i < ss_pos.size(); i++)
		in >> ss_pos[i];
	}
	
	// read in the number of defined surface cells for each position
	if (keyword == "num_defined_surcells:")
	{
	    // make sure number of user-/host-defined surcells has been sized
	    Insist (num_defined_surcells.size() != 0,
		    "Number of surface sources = 0!");

	    // input the number of user-/host-defined surcells.
	    // input file requires any preceding zeros.
	    for (int i = 0; i < num_defined_surcells.size(); i++)
		in >> num_defined_surcells[i];

	    // resize defined_surcells
	    for (int i = 0; i < num_defined_surcells.size(); i++)
		defined_surcells[i].resize(num_defined_surcells[i]);
	}

	// read in the defined surface cells list; keyword repeated for each
	// surface source with defined surface cells
	if (keyword == "defined_surcells:")
	{
	    // make sure user-/host-defined surcells has been sized
	    Insist (defined_surcells.size() != 0, 
		    "Number of surface sources = 0!");

	    // surface source index (first entry after keyword)
	    in >> data;
	    Insist (data <= defined_surcells.size(), 
		    "Defined surface source out of range!");

            Insist (defined_surcells[data-1].size() ==
		    num_defined_surcells[data-1], 
		    "Range error in defined surface source cells!");

	    // input user-/host-defined surcells
	    for (int i = 0; i < num_defined_surcells[data-1]; i++)
		in >> defined_surcells[data-1][i];
	}
    }
}

//---------------------------------------------------------------------------//
// PRIVATE MESH BUILDING IMPLEMENTATION
//---------------------------------------------------------------------------//

OS_Builder::SP_Mesh OS_Builder::build_2DMesh(SP_Coord_sys coord,  
					     Layout &layout)
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
    vf_double vertex(dimension);
    vf_int cell_pair(num_cells);

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
    SP_Mesh mesh_return(new OS_Mesh(coord, layout, vertex, cell_pair));

    // return mesh to builder
    return mesh_return;
}

//---------------------------------------------------------------------------//

OS_Builder::SP_Mesh OS_Builder::build_3DMesh(SP_Coord_sys coord, 
					     Layout &layout)
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
    vf_double vertex(dimension);
    vf_int cell_pair(num_cells);

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
    SP_Mesh mesh_return(new OS_Mesh(coord, layout, vertex, cell_pair));

    // return mesh to builder
    return mesh_return;
}

//---------------------------------------------------------------------------//
// PRIVATE COORD_SYS BUILDER IMPLEMENTATION
//---------------------------------------------------------------------------//

OS_Builder::SP_Coord_sys OS_Builder::build_Coord()
{
    // build coordinate system
    SP_Coord_sys coord;
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
    else if (coord_system == "rz" || coord_system == "RZ")
    {
	SP<XYZCoord_sys> xyzcoord(new XYZCoord_sys);
	coord = xyzcoord;
    }

    // return base class SP to a derived Coord_sys
    return coord;
}

//---------------------------------------------------------------------------//
// PRIVATE LAYOUT BUILDER IMPLEMENTATION
//---------------------------------------------------------------------------//

OS_Builder::SP_Layout OS_Builder::build_Layout(const Coord_sys &coord)
{
    // set size of new Layout
    int size = 1;
    for (int d = 0; d < coord.get_dim(); d++)
	size *= fine_edge[d].size() - 1;
    SP_Layout layout(new Layout(size));

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

//---------------------------------------------------------------------------//
// CELL ZONING FUNCTIONS (PRIVATE IMPLEMENTATION)
//---------------------------------------------------------------------------//
// Map mesh zone indices to cells

void OS_Builder::zone_mapper()
{
    // dimension of problem
    int dim = fine_edge.size();

    // determine number of zones and accumulated fine_cell data
    int num_cells = 1;
    int num_zones = 1;
    for (int d = 0; d < dim; d++)
    {
	// calculate number of zones
	num_cells *= fine_edge[d].size() - 1;
	num_zones *= coarse_edge[d].size() - 1;

	// calculate accumulated fine_cell data
        accum_cells[d].resize(fine_cells[d].size()+1);
	accum_cells[d][0] = 0;
	for (int i = 0; i < fine_cells[d].size(); i++)
	    accum_cells[d][i+1] = accum_cells[d][i] + fine_cells[d][i];
    }

    // size zone arrays
    zone.resize(num_cells);
    cell_zone.resize(num_zones);
  
    // determine cell-zone map, ie. map a fine-cell to a zone
    if (dim == 2)
        for (int j = 1; j <= fine_cells[1].size(); j++)
            for (int i = 1; i <= fine_cells[0].size(); i++)
                cell_zoner(i, j);
    else if (dim == 3)
        for (int k = 1; k <= fine_cells[2].size(); k++)
            for (int j = 1; j <= fine_cells[1].size(); j++)
                for (int i = 1; i <= fine_cells[0].size(); i++)
                    cell_zoner(i, j, k);
}
  
//---------------------------------------------------------------------------//
// Map cells to zones in 2D mesh.

void OS_Builder::cell_zoner(int iz, int jz)
{
    // match a fine-cell to a zone for 2D meshes

    // descriptive variables
    int num_xzones = coarse_edge[0].size() - 1;
    int zone_index = 1 + (iz-1) + num_xzones * (jz-1);
    int num_xcells = fine_edge[0].size() - 1;

    // loop boundaries
    int starti = 1 + accum_cells[0][iz-1];
    int endi   = accum_cells[0][iz-1] + fine_cells[0][iz-1];
    int startj = 1 + accum_cells[1][jz-1];
    int endj   = accum_cells[1][jz-1] + fine_cells[1][jz-1];

    // loop over zone and assign zone_index to those fine-cells residing in
    // the zone
    for (int j = startj; j <= endj; j++)
        for (int i = starti; i <= endi; i++)
        {
            int cell     = 1 + (i-1) + num_xcells * (j-1);
            zone[cell-1] = zone_index;
	    cell_zone[zone_index-1].push_back(cell);
        }
}

//---------------------------------------------------------------------------//
// Map cells to zones in 3D mesh.

void OS_Builder::cell_zoner(int iz, int jz, int kz)
{
    // match a fine-cell to a zone for 3D meshes

    // descriptive variables
    int num_xzones = coarse_edge[0].size() - 1;
    int num_yzones = coarse_edge[1].size() - 1;
    int zone_index = 1 + (iz-1) + num_xzones * (jz-1) + num_xzones *
	num_yzones * (kz-1);
    int num_xcells = fine_edge[0].size() - 1;
    int num_ycells = fine_edge[1].size() - 1;

    // loop boundaries
    int starti = 1 + accum_cells[0][iz-1];
    int endi   = accum_cells[0][iz-1] + fine_cells[0][iz-1];
    int startj = 1 + accum_cells[1][jz-1];
    int endj   = accum_cells[1][jz-1] + fine_cells[1][jz-1];
    int startk = 1 + accum_cells[2][kz-1];
    int endk   = accum_cells[2][kz-1] + fine_cells[2][kz-1];

    // loop over zone and assign zone_index to those fine-cells residing in
    // the zone
    for (int k = startk; k <= endk; k++)
	for (int j = startj; j <= endj; j++)
	    for (int i = starti; i <= endi; i++)
	    {
		int cell = 1 + (i-1) + num_xcells * (j-1) + num_xcells *
		    num_ycells * (k-1);
		zone[cell-1] = zone_index;
		cell_zone[zone_index-1].push_back(cell);
	    }
}

//---------------------------------------------------------------------------//
// CALCULATE DEFINED SURFACE CELLS
//---------------------------------------------------------------------------//
/*!  
 * \brief Use built global mesh to calculate a list of defined surface
 * cells.
 *
 * This function uses the global mesh to convert a global boundary surface
 * source into a list of explicit, global cells along the boundary.  For
 * example, a "lox" boundary could be defined as a surface source.  This
 * description will be converted into a list of all global cells along the
 * low-x boundary.
 *
 * Surface source cell lists are only calculated for surface sources that do
 * not already have user-defined lists of cells.
 */
void OS_Builder::calc_defined_surcells()
{
    Require (mesh);
    Require (mesh->full_Mesh());

    // define number of surface sources
    int num_ss = ss_pos.size();

    // exit if no surface sources
    if (num_ss == 0)
	return;

    // check the defined surface source list
    Check (defined_surcells.size() == num_ss);

    // check to see if we need to do any work
    for (int i = 0; i < num_ss; i++)
    {
	if (defined_surcells[i].size() == 0)
	{
	    // get a vector of surface cells along the boundary
	    defined_surcells[i] = mesh->get_surcells(ss_pos[i]);
	    Check (defined_surcells[i].size() > 0);
	}
	else if (defined_surcells[i].size() > 0)
	{
	    Check (mesh->check_defined_surcells
		   (ss_pos[i], defined_surcells[i]));
	    continue;
	}
	else
	{
	    Insist (0, "Number of defined surface source cells < 0!");
	}
    }
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of OS_Builder.cc
//---------------------------------------------------------------------------//
