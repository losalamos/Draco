//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Interface.cc
// B.T. Adams (bta@lanl.gov)
// 18 May 99
//---------------------------------------------------------------------------//
// @> CAR_CU_Interface class implementation file (developed from 
// OS_Interface.hh)
//---------------------------------------------------------------------------//

#include "CAR_CU_Interface.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <map>
#include <iterator>

namespace rtt_imc 
{

using std::fill;
using std::cout;
using std::endl;
using std::sort;
using std::multimap;
using std::make_pair;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public Parser member functions
//---------------------------------------------------------------------------//

SP<RTT_Format> CAR_CU_Interface::parser()
{
  // open input file, ifstream object requires C-style string
    const char *file = input_file.c_str();
    ifstream input(file);
    SP<RTT_Format> spRTT_Format;

  // make sure file exists
    Insist (input, "You must supply an input file!");

  // call CAR_CU_Mesh Parser
    spRTT_Format = parser_Mesh(input);

  // call Opacities Parser
    parser_Opacity(input);

  // call Source Parser
    parser_Source(input);

    return spRTT_Format;
}

//---------------------------------------------------------------------------//
// private Parser member functions for CAR_CU_Mesh build
//---------------------------------------------------------------------------//
// parse the Mesh

SP<RTT_Format> CAR_CU_Interface::parser_Mesh(ifstream &in)
{
  // CAR_CU_Mesh Parser

  // determine coord_sys
    string keyword;
    while (keyword != "end-title")
    {
      // test that we have not reached end-of-file
	Insist (!in.eof(), "You did not supply an end-title!");

      // do title block input
	in >> keyword;
	if (keyword == "coord:")
	{
	    in >> coord_system;
	    if (coord_system == "xyz")
	        coord_system = "XYZ";
	    else if (coord_system != "XYZ")
	        Insist (0, "Only XYZ geometry currently supported!")
	}
	if (keyword == "rtt_file:")
	    in >> rtt_mesh_file;
	if (keyword == "surface_file:")
	    in >> surface_file ;
    }
    Insist(rtt_mesh_file !="undefined",
	       "You did not supply a RTT_Mesh format input file!");

    // Set a flag to renumber the nodes, sides, and cells into ascending
    // order based on coordinates with x spinning fastest, followed by y,
    // and finally z. The mapping of faces has been changed from milagro's 
    // OS_Mesh definition (-z = 4, -y = 2, -x = 0, x = 1, y = 3, z = 5) to the
    // rttMesh definition (-z = 0, -y = 1, -x = 2, x = 3, y = 4, z = 5). Read 
    // the data from the rtt-format file.
    bool renumber = true;
    SP<RTT_Format> rttMesh(new RTT_Format(rtt_mesh_file, renumber));

    // Set up some useful RTT Format variables
    // total number of nodes.
    int ncells = rttMesh->get_dims_ncells();
    // problem dimension (2D or 3D).

    // Assign a material zone to each cell.
    zone.resize(ncells);
    int matl_flag = rttMesh->get_cell_flags_material_flag_number();
    for (int cell_number = 0; cell_number < ncells; cell_number++)
        zone[cell_number] = rttMesh->get_cells_flags(cell_number,matl_flag);

    // resize the material zone vector according to the number of materials
    mat_zone.resize(rttMesh->get_dims_ncell_flags(matl_flag));

    return rttMesh;
}

//---------------------------------------------------------------------------//
// parse and assign the opacity and mat_state variables

void CAR_CU_Interface::parser_Opacity(ifstream &in)
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
	    Insist(mat_zone.size() == data, 
		   "Number of materials and zones must be equal!");
	    density.resize(data);
	    kappa.resize(data);
	    kappa_thomson.resize(data);
	    temperature.resize(data);
	    specific_heat.resize(data);
	}
	if (keyword == "mat:")
	{
	  // Check that material data arrays have been sized previously
	    Check (density.size() != 0);
	    Check (density.size() == kappa.size());
	    Check (density.size() == kappa_thomson.size());
	    Check (density.size() == temperature.size());
	    Check (density.size() == specific_heat.size());
	  // input the material arrays
	    for (int i = 0; i < density.size(); i++)
	    {
		in >> data;
		in >> density[data-1];
		in >> kappa[data-1];
		in >> kappa_thomson[data-1];
		in >> temperature[data-1];
		in >> specific_heat[data-1];
	    }
	}
	if (keyword == "implicitness:")
	    in >> implicitness;
	if (keyword == "analytic_opacity:")
	    in >> analytic_opacity;
        if (keyword == "analytic_sp_heat:")
            in >> analytic_sp_heat;
    }    

  // make sure we have gotten a fleck factor
    Insist(implicitness >= 0 && implicitness <= 1, 
	   "You must specify a proper Implicitness factor!");	
    Insist(density.size() > 0, "You must specify at least 1 Mat.!");
}


//---------------------------------------------------------------------------//
// private Parser member functions for Source_Init<MT> build
//---------------------------------------------------------------------------//
// parse input for Source_Init

void CAR_CU_Interface::parser_Source(ifstream &in)
{
  // define number of zones and resize Source data arrays
    int num_zones = mat_zone.size();
    evol_ext.resize(num_zones);
    rad_temp.resize(num_zones);
    rad_source.resize(num_zones);
    fill(evol_ext.begin(), evol_ext.end(), 0.0);
    fill(rad_temp.begin(), rad_temp.end(), 0.0);
    fill(rad_source.begin(), rad_source.end(), 0.0);

  // parse the source
    zone_source_parser(in);
}

//---------------------------------------------------------------------------//
// parse and assign the source variables
    
void CAR_CU_Interface::zone_source_parser(ifstream &in)
{
  // Source Parser

  // input keywords
    string keyword;
    int    data;
    bool   have_rad_source = false;

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
	if (keyword == "rad_source:")
	{
	    for (int i = 0; i < rad_source.size(); i++)
	    {
		in >> rad_source[i];
		if (rad_source[i] > 0.0)
		    have_rad_source = true;
	    }
	}
	if (keyword == "rad_s_tend:")
	    in >> rad_s_tend;
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
	if (keyword == "print_frequency:")
	    in >> print_f;
	if (keyword == "buffer_size:")
	    in >> buffer;
	if (keyword == "seed:")
	    in >> seed;
    }   
  
  // do some assertions on the source variables
    Insist (delta_t > 0, "The timestep must be > 0!");
    Insist (npmax > 0, "The npmax must be > 0!");
    Insist (npnom > 0, "The npnom must be > 0!");
    Insist (capacity > 0, "The capacity must be > 0!");
    Insist (max_cycle > 0, "The max_cycle must be > 0!");
    Insist (buffer > 0, "The buffer must be > 0!");
    if (have_rad_source)
	Insist (rad_s_tend > 0, "You have a radiation source of zero time!");
    if (rad_s_tend > 0)
	Insist (have_rad_source, 
		"Duration defined for zero radiation source!");
}

//---------------------------------------------------------------------------//
// Zone - to - cell mapping accessor functions
//---------------------------------------------------------------------------//
// map kappa to cell based arrays

vector<double> CAR_CU_Interface::get_kappa() const
{
  // make a return vector of the proper size
    vector<double> cell_k(zone.size());

  // assign the values to cell_k based on the zone and material
    for (int cell = 1; cell <= cell_k.size(); cell++)
	cell_k[cell-1] = kappa[mat_zone[zone[cell-1]-1]-1];

  // return cell_k
    return cell_k;
}

//---------------------------------------------------------------------------//
// map kappa_thomson to cell-based arrays

vector<double> CAR_CU_Interface::get_kappa_thomson() const
{
  // make a return vector of the proper size
    vector<double> cell_kt(zone.size());

  // assign the values to cell_kt based on the zone and material
    for (int cell = 1; cell <= cell_kt.size(); cell++)
	cell_kt[cell-1] = kappa_thomson[mat_zone[zone[cell-1]-1]-1];

  // return cell_k
    return cell_kt;
}

//---------------------------------------------------------------------------//
// map density to cell based arrays

vector<double> CAR_CU_Interface::get_density() const
{
  // make a return vector of the proper size
    vector<double> cell_den(zone.size());

  // assign the values to cell_den based on the zone and material
    for (int cell = 1; cell <= cell_den.size(); cell++)
	cell_den[cell-1] = density[mat_zone[zone[cell-1]-1]-1];

  // return cell_den
    return cell_den;
}

//---------------------------------------------------------------------------//
// map temperature to cell based arrays

vector<double> CAR_CU_Interface::get_temperature() const
{
  // make a return vector of the proper size
    vector<double> cell_t(zone.size());

  // assign the values to cell_t based on the zone and material
    for (int cell = 1; cell <= cell_t.size(); cell++)
	cell_t[cell-1] = temperature[mat_zone[zone[cell-1]-1]-1];

  // return cell_t
    return cell_t;
}

//---------------------------------------------------------------------------//
// map specific heat to cell based arrays

vector<double> CAR_CU_Interface::get_specific_heat() const
{
  // make a return vector of the proper size
    vector<double> cell_cv(zone.size());

  // assign the values to cell_cv based on the zone and material
    for (int cell = 1; cell <= cell_cv.size(); cell++)
	cell_cv[cell-1] = specific_heat[mat_zone[zone[cell-1]-1]-1];

  // return cell_cv
    return cell_cv;
}

//---------------------------------------------------------------------------//
// map volume source to cell based arrays

vector<double> CAR_CU_Interface::get_evol_ext() const
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
// map radiation source to cell based arrays

vector<double> CAR_CU_Interface::get_rad_source() const
{
  // make a return vector of the proper size
    vector<double> cell_rsrc(zone.size());

  // assign cell values of rad source based on zonal values
    for (int cell = 1; cell <= cell_rsrc.size(); cell++)
	cell_rsrc[cell-1] = rad_source[zone[cell-1]-1];

  // return cell_evol
    return cell_rsrc;
}

//---------------------------------------------------------------------------//
// map rad temp to cell based arrays

vector<double> CAR_CU_Interface::get_rad_temp() const
{
  // make a return vector of the proper size
    vector<double> cell_rad(zone.size());

  // assign the values to cell_rad based on the zone
    for (int cell = 1; cell <= cell_rad.size(); cell++)
	cell_rad[cell-1] = rad_temp[zone[cell-1]-1];

  // return rad_temp
    return cell_rad;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of CAR_CU_Interface.cc
//---------------------------------------------------------------------------//
