//----------------------------------*-C++-*----------------------------------//
// OS_Interface.cc
// Thomas M. Evans
// Mon Feb 23 17:22:21 1998
//---------------------------------------------------------------------------//
// @> OS_Interface class implementation file
//---------------------------------------------------------------------------//

#include "imc/OS_Interface.hh"
#include "ds++/Assert.hh"
#include <algorithm>

IMCSPACE

using std::fill;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public Parser member functions
//---------------------------------------------------------------------------//

void OS_Interface::parser()
{
  // open input file, ifstream object requires C-style string
    const char *file = input_file.c_str();
    ifstream input(file);

  // make sure file exists
    Insist (input, "You must supply an input file!");

  // call OS_Mesh Parser
    parser_Mesh(input);

  // call Opacities Parser
    parser_Opacity(input);

  // call Source Parser
    parser_Source(input);
}

//---------------------------------------------------------------------------//
// private Parser member functions for OS_Mesh build
//---------------------------------------------------------------------------//
// parse the Mesh

void OS_Interface::parser_Mesh(ifstream &in)
{
  // OS_Mesh Parser

  // determine coord_sys
    string keyword;
    while (keyword != "end-title")
    {
      // test that we have not reached end-of-file
	Insist (!in.eof(), "You did not supply an end-title!");

      // do title block input
	in >> keyword;
	if (keyword == "coord:")
	    in >> coord_system;
    }

  // check to make sure an appropriate coord_sys is called
    Check (coord_system == "xy" || coord_system == "XY" || 
	   coord_system == "rz" || coord_system == "RZ" ||
	   coord_system == "XYZ" || coord_system == "xyz");

  // call appropriate sub-parser
    if (coord_system == "xy" || coord_system == "XY" || coord_system == "rz"
	|| coord_system == "RZ")
    {
	fine_cells.resize(2);
	accum_cells.resize(2);
	coarse_edge.resize(2);
	fine_edge.resize(2);
	bnd_cond.resize(4);
	parser2D(in);
    }
    else if (coord_system == "xyz" || coord_system == "XYZ")
    {
	fine_cells.resize(3);
	accum_cells.resize(3);
	coarse_edge.resize(3);
	fine_edge.resize(3);
	bnd_cond.resize(6);
	parser3D(in);
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

//---------------------------------------------------------------------------//
// parse a 2D mesh

void OS_Interface::parser2D(ifstream &in)
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
	Insist (!in.eof(), "You did not specify a proper end-mesh!");

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

//---------------------------------------------------------------------------//
// parse a 3D mesh

void OS_Interface::parser3D(ifstream &in)
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
	Insist (!in.eof(), "You did not specify a proper end-mesh!");
     
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
// private Parser member functions for Opacity<MT> build
//---------------------------------------------------------------------------//
// parse input for Opacity and Mat_State class building

void OS_Interface::parser_Opacity(ifstream &in)
{ 
  // determine size of mesh and initialize zone and zone_map arrays
    int dim       = fine_edge.size();
    int num_cells = 1;
    int num_zones = 1;
    for (int d = 0; d < dim; d++)
    {
        num_cells *= fine_edge[d].size() - 1;
        num_zones *= coarse_edge[d].size() - 1;
    }
    zone.resize(num_cells);
    mat_zone.resize(num_zones);

  // map fine cells to zones
    zone_mapper();

  // read zone map
    zone_opacity_parser(in);
}

//---------------------------------------------------------------------------//
// parse and assign the opacity and mat_state variables

void OS_Interface::zone_opacity_parser(ifstream &in)
{
  // Opacity Parser

  // input keywords
    string keyword;
    int    data;

  // determine zone map
    while (keyword != "end-mat")
    {
      // test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-mat!");
	
      // do input
	in >> keyword;
	if (keyword == "zonemap:")
	    for (int i = 0; i < mat_zone.size(); i++)
	    {
		in >> mat_zone[i];
		Insist(mat_zone[i] >= 0, "Materials in zonemap <= 0!");
	    }
	if (keyword == "num_materials:")
	{
	    in >> data;
	    density.resize(data);
	    kappa.resize(data);
	    temperature.resize(data);
	    specific_heat.resize(data);
	}
	if (keyword == "mat:")
	{
	  // Check that material data arrays have been sized previously
	    Check (density.size() != 0);
	    Check (density.size() == kappa.size());
	    Check (density.size() == temperature.size());
	    Check (density.size() == specific_heat.size());
	  // input the material arrays
	    for (int i = 0; i < density.size(); i++)
	    {
		in >> data;
		in >> density[data-1];
		in >> kappa[data-1];
		in >> temperature[data-1];
		in >> specific_heat[data-1];
	    }
	}
	if (keyword == "implicitness:")
	    in >> implicitness;
	if (keyword == "analytic_opacity:")
	    in >> analytic_opacity;
    }    

  // make sure we have gotten a fleck factor
    Insist(implicitness >= 0 && implicitness <= 1, 
	   "You must specify a proper Implicitness factor!");	
    Insist(density.size() > 0, "You must specify at least 1 Mat.!");
}

//---------------------------------------------------------------------------//
// map fine cells to zones

void OS_Interface::zone_mapper()
{
  // dimension of mesh
    int dim = fine_edge.size();

  // calculate accumulated fine_cell array
    for (int d = 0; d < dim; d++)
    {
        accum_cells[d].resize(fine_cells[d].size()+1);
	accum_cells[d][0] = 0;
	for (int i = 0; i < fine_cells[d].size(); i++)
	    accum_cells[d][i+1] = accum_cells[d][i] + fine_cells[d][i];
    }
  
  // determine cell-zone map, ie. map a fine-cell to a zone
    if (dim == 2)
        for (int j = 1; j <= fine_cells[1].size(); j++)
            for (int i = 1; i <= fine_cells[0].size(); i++)
                cell_zone(i, j);
    else if (dim == 3)
        for (int k = 1; k <= fine_cells[2].size(); k++)
            for (int j = 1; j <= fine_cells[1].size(); j++)
                for (int i = 1; i <= fine_cells[0].size(); i++)
                    cell_zone(i, j, k);
}
  
//---------------------------------------------------------------------------//

void OS_Interface::cell_zone(int iz, int jz)
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
        }
}

//---------------------------------------------------------------------------//

void OS_Interface::cell_zone(int iz, int jz, int kz)
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
	    }
}

//---------------------------------------------------------------------------//
// private Parser member functions for Source_Init<MT> build
//---------------------------------------------------------------------------//
// parse input for Source_Init

void OS_Interface::parser_Source(ifstream &in)
{
  // define number of zones and resize Source data arrays
    int num_zones = mat_zone.size();
    evol_ext.resize(num_zones);
    rad_temp.resize(num_zones);
    fill(evol_ext.begin(), evol_ext.end(), 0.0);
    fill(rad_temp.begin(), rad_temp.end(), 0.0);

  // parse the source
    zone_source_parser(in);
}

//---------------------------------------------------------------------------//
// parse and assign the source variables
    
void OS_Interface::zone_source_parser(ifstream &in)
{
  // Source Parser

  // input keywords
    string keyword;
    int    data;

  // determine zone map
    while (keyword != "end-source")
    {
      // test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-source!");
	
      // do input
	in >> keyword;
	if (keyword == "vol_source:")
	    for (int i = 0; i < evol_ext.size(); i++)
		in >> evol_ext[i];
	if (keyword == "num_ss:")
	{
	    in >> data;
	    ss_pos.resize(data);
	    ss_temp.resize(data);
	    fill(ss_temp.begin(), ss_temp.end(), 0.0);
	}
	if (keyword == "sur_source:")
	{
	  // make sure that surface source has been sized
	    Check (ss_pos.size() != 0);
	  // input the surface source position
	    for (int i = 0; i < ss_pos.size(); i++)
		in >> ss_pos[i];
	}
	if (keyword == "sur_temp:")
	{
	  // make sure that surface source has been sized
	    Check (ss_temp.size() != 0);
	  // input the surface source temps
	    for (int i = 0; i < ss_temp.size(); i++)
		in >> ss_temp[i];
	}
	if (keyword == "ss_dist:")
	{
	  // make sure a ss has been defined
	    Check (ss_pos.size() != 0);
	    in >> ss_dist;
	}
	if (keyword == "rad_temp:")
	    for (int i = 0; i < rad_temp.size(); i++)
		in >> rad_temp[i];
	if (keyword == "timestep:")
	    in >> delta_t;
	if (keyword == "npmax:")
	    in >> npmax;
	if (keyword == "npnom:")
	    in >> npnom;
	if (keyword == "dnpdt:")
	    in >> dnpdt;
	if (keyword == "capacity:")
	    in >> capacity;
	if (keyword == "max_cycle:")
	    in >> max_cycle;
    }   
  
  // do some assertions on the source variables
    Insist (delta_t > 0, "The timestep must be > 0!");
    Insist (npmax > 0, "The npmax must be > 0!");
    Insist (npnom > 0, "The npnom must be > 0!");
    Insist (capacity > 0, "The capacity must be > 0!");
    Insist (max_cycle > 0, "The max_cycle must be > 0!");
}

//---------------------------------------------------------------------------//
// map volume source to cell based arrays

vector<double> OS_Interface::get_evol_ext() const
{
  // make a return vector of the proper size
    vector<double> cell_evol(zone.size());

  // assign the values to cell_evol based on the zone
    for (int cell = 1; cell <= cell_evol.size(); cell++)
	cell_evol[cell-1] = evol_ext[zone[cell-1]-1];

  // return cell_evol
    return cell_evol;
}

//---------------------------------------------------------------------------//
// map rad temp to cell based arrays

vector<double> OS_Interface::get_rad_temp() const
{
  // make a return vector of the proper size
    vector<double> cell_rad(zone.size());

  // assign the values to cell_rad based on the zone
    for (int cell = 1; cell <= cell_rad.size(); cell++)
	cell_rad[cell-1] = rad_temp[zone[cell-1]-1];

  // return rad_temp
    return cell_rad;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of OS_Interface.cc
//---------------------------------------------------------------------------//
