//----------------------------------*-C++-*----------------------------------//
// OS_Parser.cc
// Thomas M. Evans
// Mon Feb 23 17:22:21 1998
//---------------------------------------------------------------------------//
// @> OS_Parser class implementation file
//---------------------------------------------------------------------------//

#include "OS_Parser.hh"
#include <cassert>

IMCSPACE

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public Parser member functions
//---------------------------------------------------------------------------//
void OS_Parser::Parser()
{
  // open input file, ifstream object requires C-style string
    const char *file = input_file.c_str();
    ifstream input(file);

  // make sure file exists
    assert (input);

  // call OS_Mesh Parser
    Parser_mesh(input);

  // call Opacities Parser
    Parser_opacity(input);

  // call Source Parser
}

//---------------------------------------------------------------------------//
// private Parser member functions for OS_Mesh build
//---------------------------------------------------------------------------//
void OS_Parser::Parser_mesh(ifstream &in)
{
  // OS_Mesh Parser

  // determine coord_sys
    string keyword;
    while (keyword != "end-title")
    {
      // test that we have not reached end-of-file
	assert (!in.eof());

      // do title block input
	in >> keyword;
	if (keyword == "coord:")
	    in >> coord_system;
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
	accum_cells.resize(2);
	coarse_edge.resize(2);
	fine_edge.resize(2);
	bnd_cond.resize(4);
	Parser2D(in);
    }
    else if (coord_system == "xyz" || coord_system == "XYZ")
    {
	fine_cells.resize(3);
	accum_cells.resize(3);
	coarse_edge.resize(3);
	fine_edge.resize(3);
	bnd_cond.resize(6);
	Parser3D(in);
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
void OS_Parser::Parser2D(ifstream &in)
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

//---------------------------------------------------------------------------//
void OS_Parser::Parser3D(ifstream &in)
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
// private Parser member functions for Opacity<MT> build
//---------------------------------------------------------------------------//
void OS_Parser::Parser_opacity(ifstream &in)
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
    Zone_mapper();
}

//---------------------------------------------------------------------------//
void OS_Parser::Zone_mapper()
{
  // dimension of mesh
    int dim = fine_edge.size();

  // calculate accumulated fine_cell array
    for (int d = 0; d < dim; d++)
    {
	accum_cells[d][0] = 0;
	for (int i = 0; i < fine_cells[d].size(); i++)
	    accum_cells[d][i+1] = accum_cells[d][i] + fine_cells[d][i];
    }
  
  // determine cell-zone map, ie. map a fine-cell to a zone
    if (dim == 2)
        for (int j = 1; j <= fine_cells[1].size(); j++)
            for (int i = 1; i <= fine_cells[0].size(); i++)
                Cell_zone(i, j);
    else if (dim == 3)
        for (int k = 1; k <= fine_cells[2].size(); k++)
            for (int j = 1; j <= fine_cells[1].size(); j++)
                for (int i = 1; i <= fine_cells[0].size(); i++)
                    Cell_zone(i, j, k);
}
  
//---------------------------------------------------------------------------//
void OS_Parser::Cell_zone(int iz, int jz)
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
void OS_Parser::Cell_zone(int iz, int jz, int kz)
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

CSPACE

//---------------------------------------------------------------------------//
//                              end of OS_Parser.cc
//---------------------------------------------------------------------------//
