//----------------------------------*-C++-*----------------------------------//
// Shadow_CAR_CU_Interface.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
//---------------------------------------------------------------------------//
// @> Shadow_CAR_CU_Interface interface file
//---------------------------------------------------------------------------//

#ifndef __mc_Shadow_CAR_CU_Interface_cc__
#define __mc_Shadow_CAR_CU_Interface_cc__

#include "CAR_CU_Interface.hh"
#include "RTT_Format.hh"
#include "Shadow_Opaque_Pointers.hh"
#include <iostream>

//===========================================================================//
// Shadow_CAR_CU_Interface - 
//
// Purpose : Provides shadow interface functions to the Continuous Adaptive 
// Refinement Cartesion Unstructured Mesh Interface Class for use with 
// Fortran 90 codes.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_imc 
{
using std::cout;
using std::endl;

// draco components
using dsxx::SP;
using rtt_format::RTT_Format;

extern "C" 
{
//---------------------------------------------------------------------------//
// CAR_CU_Interface F90 to C++ flat interface functions
//---------------------------------------------------------------------------//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
//
    // Construct a CAR_CU_Interface class from a Fortran 90 program call. This
    // also constructs a RTT_Format class object and parses both the input
    // deck and the RTT Format mesh file specified therein. The addresses of
    // both the new CAR_CU_Interface and RTT_Format class objects are set.
    void construct_car_cu_interface_(long & self, char * file, int & verbosity,
				     long & rttFormat)
    {
        string infile = file;
	bool verbose = verbosity;

	// Construct a new CAR_CU_Interface class object.
	SP<CAR_CU_Interface> interface(new CAR_CU_Interface(infile, verbose));

	// Construct a new RTT_Format class object and parse the input files.
	SP<RTT_Format> rttMesh = interface->parser();
	if (verbose)
	    cout << " ** Read input file " << infile << endl;

	// return the addresses of the new CAR_CU_Interface (self) and 
	// RTT_Format class (rttFormat) objects.
	self = opaque_pointers<CAR_CU_Interface>::insert(interface);
	rttFormat = opaque_pointers<RTT_Format>::insert(rttMesh);

    }

    // Destroy a CAR_CU_Interface class object from a Fortran 90 program call.
    void destruct_car_cu_interface_(long & self)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// destroy the CAR_CU_Interface class object by assigning this SP to 
        // a null SP
	interface = SP<CAR_CU_Interface>();
	Ensure (!interface);

	// remove the opaque pointer to the CAR_CU_Interface class object.
	opaque_pointers<CAR_CU_Interface>::erase(self);
    }

//===========================================================================//
// General CAR_CU_Interface accessor functions
//===========================================================================//
    // Return the number of sets of grouped surface source cells.
    void get_car_cu_ss_pos_size_(long & self, long & ss_pos_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	int isize = interface->get_ss_pos_size();

	// Cast the int variable to long
	ss_pos_size = static_cast<long>(isize);
    }

    // Return the number of grouped surface source cells in a given set.
    void get_car_cu_ss_cells_size_(long & self, long & surface, 
				   long & ss_cells_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isurface = static_cast<int>(surface);
	Insist(isurface > 0 && isurface <= interface->get_ss_pos_size(), 
	    "Invalid surface number passed to get_car_cu_ss_cells_size_!");

	int isize = interface->get_ss_cells_size(isurface);

	// Cast the int variable to long
	ss_cells_size = static_cast<long>(isize);
    }

    // Return the position (lox, hix, etc.) of all of the grouped surface 
    // source cell sets.
    void get_car_cu_ss_pos_(long & self, char * pos, long & ss_pos_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isize = static_cast<int>(ss_pos_size);

	Insist(isize == 3 * interface->get_ss_pos_size(), 
	    "Invalid ss position size passed to get_car_cu_ss_pos_!");

	vector<string> ss_pos = interface->get_ss_pos();

	for (int surf = 0; surf < ss_pos.size(); surf++)
	{
	    for (int ind = 0; ind < ss_pos[surf].size(); ind ++)
	    {
	        * pos = ss_pos[surf][ind];
		++pos;
	    }
	}
    }

    // Return the position (lox, hix, etc.) of a set of grouped surface source 
    // cells.
    void get_car_cu_ss_cells_pos_(long & self, long & surface, char * pos)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isurface = static_cast<int>(surface);

	Insist(isurface > 0 && isurface <= interface->get_ss_pos_size(), 
	    "Invalid surface number passed to get_car_cu_ss_cells_pos_!");

	string ss_pos = interface->get_ss_pos(isurface);

	for (int ind = 0; ind < ss_pos.size(); ind ++)
	{
	    * pos = ss_pos[ind];
	    ++pos;
	}
    }

    // Return the surface source angular distribution.
    void get_car_cu_ss_dist_(long & self, char * distribution)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	string ss_dist = interface->get_ss_dist();

	for (int ind = 0; ind < ss_dist.size(); ind ++)
	{
	    * distribution = ss_dist[ind];
	    ++distribution;
	}
    }
    // Return the temperature of all of the grouped surface source cell sets.
    void get_car_cu_ss_temp_(long & self, double & temp, long & ss_temp_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isize = static_cast<int>(ss_temp_size);
	double * data_array = & temp;

	Insist(isize == interface->get_ss_pos_size(), 
	    "Invalid ss temperature size passed to get_car_cu_ss_temp_!");

	vector<double> ss_temp = interface->get_ss_temp();

	for (int surf = 0; surf < ss_temp.size(); surf++)
	{
	        * data_array = ss_temp[surf];
		++data_array;
	}
    }

    // Return the temperature of a set of grouped surface source cells.
    void get_car_cu_ss_cells_temp_(long & self, long & surface, double & temp)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isurface = static_cast<int>(surface);

	Insist(isurface > 0 && isurface <= interface->get_ss_pos_size(), 
	    "Invalid surface number passed to get_car_cu_ss_cells_temp_!");

	temp = interface->get_ss_temp(isurface);
    }

    // Return all of the defined surface source cell sets.
    void get_car_cu_ss_cells_(long & self, long & ss_set, long & ss_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variables to int
	int isize = static_cast<int>(ss_size);
	long * data_array = & ss_set;

	int size_check = 0;
	for (int surf = 1; surf <= interface->get_ss_pos_size(); surf++)
	    size_check += interface->get_ss_cells_size(surf);

	Insist(isize == size_check, 
	    "Invalid surface source cells size passed to get_ss_cells_!");

	vector<vector<int> > iss_set = interface->get_defined_surcells();

	for (int surf = 0; surf < iss_set.size(); surf++)
	{
	    for (int cell = 0; cell < iss_set[surf].size(); cell++)
	    {
	        * data_array = static_cast<long>(iss_set[surf][cell]);
		++data_array;
	    }
	}
    }

    // Return the defined surface source cells in a given set.
    void get_car_cu_ss_cells_set_(long & self, long & surface, long & ss_set, 
				  long & ss_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variables to int
	int isurface = static_cast<int>(surface);
	int isize = static_cast<int>(ss_size);
	long * data_array = & ss_set;

	Insist(isurface > 0 && isurface <= interface->get_ss_pos_size(), 
	    "Invalid surface number passed to get_car_cu_ss_cell_set_!");
	Insist(isize == interface->get_ss_cells_size(isurface), 
	    "Invalid surface source cells size passed to get_ss_cell_set_!");

	vector<int> iss_set = interface->get_defined_surcells(isurface);

	for (int cell = 0; cell < iss_set.size(); cell++)
	{
	    * data_array = static_cast<long>(iss_set[cell]);
	    ++data_array;
	}
    }

    // Return the mesh volumetric sources for all cells.
    void get_car_cu_vol_src_(long & self, double & src, long & src_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isrc_size = static_cast<int>(src_size);
	double * data_array = & src;

	Insist(isrc_size == interface->get_zone_size(), 
	    "Invalid volume source size passed to get_car_cu_vol_src_!");

	vector<double> vol_src = interface->get_evol_ext();

	for (int cell = 0; cell < vol_src.size(); cell++)
	{
	        * data_array = vol_src[cell];
		++data_array;
	}
    }

    // Return the volumetric source for a single cell.
    void get_car_cu_cell_vol_src_(long & self, long & cell, double & src)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int icell = static_cast<int>(cell);

	Insist(icell > 0 && icell <= interface->get_zone_size(), 
	    "Invalid cell number passed to get_car_cu_cell_vol_src_!");

	src = interface->get_evol_ext(icell);
    }

    // Return the mesh radiation sources for all cells.
    void get_car_cu_rad_src_(long & self, double & src, long & src_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isrc_size = static_cast<int>(src_size);
	double * data_array = & src;

	Insist(isrc_size == interface->get_zone_size(), 
	    "Invalid radiation source size passed to get_car_cu_rad_src_!");

	vector<double> rad_src = interface->get_rad_source();

	for (int cell = 0; cell < rad_src.size(); cell++)
	{
	        * data_array = rad_src[cell];
		++data_array;
	}
    }

    // Return the radiation source for a single cell.
    void get_car_cu_cell_rad_src_(long & self, long & cell, double & src)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int icell = static_cast<int>(cell);

	Insist(icell > 0 && icell <= interface->get_zone_size(), 
	    "Invalid cell number passed to get_car_cu_cell_rad_src_!");

	src = interface->get_rad_source(icell);
    }

    // Return the cut-off time for a radiation source.
    void get_car_cu_rad_s_tend_(long & self, double & time)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	time = interface->get_rad_s_tend();
    }

    // Return the mesh analytic specific heat type
    void get_car_cu_analy_opacity_(long & self, char * analy_opacity)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	string opacity = interface->get_analytic_sp_heat();

	for (int index = 0; index < opacity.size(); index++)
	{
	    * analy_opacity = opacity[index];
	    ++analy_opacity;
	}
    }

    // Return the mesh analytic specific heat type
    void get_car_cu_analy_sp_heat_(long & self, char * analy_spec_heat)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	string spec_heat = interface->get_analytic_sp_heat();

	for (int index = 0; index < spec_heat.size(); index++)
	{
	    * analy_spec_heat = spec_heat[index];
	    ++analy_spec_heat;
	}
    }

    // Return the mesh implicitness factor (Fleck's alpha)
    void get_car_cu_implicit_(long & self, long & implicit)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	implicit = static_cast<long>(interface->get_implicit());

    }

    // Return the initial time step size
    void get_car_cu_time_step_(long & self, double & time_step)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	time_step = interface->get_delta_t();

    }

    // Return the processor capacity (cells/processor).
    void get_car_cu_capacity_(long & self, long & capacity)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	capacity = static_cast<long>(interface->get_capacity());

    }

    // Return the number of cycles to run.
    void get_car_cu_num_cycles_(long & self, long & num_cycles)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	num_cycles = static_cast<long>(interface->get_max_cycle());

    }

    // Return the print frequency.
    void get_car_cu_print_freq_(long & self, long & print_freq)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	print_freq = static_cast<long>(interface->get_printf());

    }

} // end extern "C"


} // end namespace rtt_imc

#endif                          // __mc_Shadow_CAR_CU_Interface_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_CAR_CU_Interface.cc
//---------------------------------------------------------------------------//
