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
    void get_car_cu_ss_pos_(long & self, long & surface, char * pos, 
			    long & ss_pos_size)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	// Cast the long variable to int
	int isurface = static_cast<int>(surface);
	int isize = static_cast<int>(ss_pos_size);

	Insist(isurface > 0 && isurface <= interface->get_ss_pos_size(), 
	    "Invalid surface number passed to get_car_cu_ss_pos_!");
	Insist(isize == interface->get_ss_pos_size(), 
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
	    "Invalid surface source size passed to get_ss_cells_!");

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
	    "Invalid surface source size passed to get_ss_cell_set_!");

	vector<int> iss_set = interface->get_defined_surcells(isurface);

	for (int cell = 0; cell < iss_set.size(); cell++)
	{
	    * data_array = static_cast<long>(iss_set[cell]);
	    ++data_array;
	}
    }

} // end extern "C"


} // end namespace rtt_imc

#endif                          // __mc_Shadow_CAR_CU_Interface_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_CAR_CU_Interface.cc
//---------------------------------------------------------------------------//
