//----------------------------------*-C++-*----------------------------------//
// Shadow_Interface.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
/*! 
 * \file   amr_mesh/Shadow_Interface.cc
 * \author B.T. Adams
 * \date   Mon 27 Sep 10:33:26 1999
 * \brief  Provides the C++ side of the shadow object interface functions to
 *         the Continuous Adaptive Refinement Cartesion Unstructured Mesh 
 *         Interface class for use with Fortran 90 codes. The complimentary 
 *         Fortran 90 shadow object interface functions that reference the 
 *         functions herein are provided in amr_mesh_fort/Shadow_Interface.f90.
 *         An example code is also provide via link to illustrate the usage
 *         of all of the shadow object interface functions to the amr_mesh 
 *         package from a Fortran 90 code.
 *
 *\sa Mark G. Gray, Randy M. Roberts, and Tom Evans, Scientific Programming,
 *   "Shadow-Object Interface Between Fortran 95 and C++", March-April 1999.
 */
//---------------------------------------------------------------------------//
// @> Shadow_Interface interface file
//---------------------------------------------------------------------------//

#ifndef __amr_Shadow_Interface_cc__
#define __amr_Shadow_Interface_cc__

#include "Interface.hh"
#include "meshReaders/RTT_Format.hh"
#include "ds++/opaquePointers.hh"
#include <iostream>

//===========================================================================//
// Shadow_Interface - 
//
// Purpose : Provides shadow object interface functions to the Continuous 
//           Adaptive Refinement Cartesion Unstructured Mesh Interface Class
//           for use with Fortran 90 codes.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_amr 
{
using std::cout;
using std::endl;

// draco components
using rtt_dsxx::SP;
using rtt_dsxx::opaque_pointers;
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
/*!
 * \brief Shadow object that constructs a CAR_CU_Interface class object from
 *        a Fortran 90 program call. This also constructs a RTT_Format class
 *        object and parses both the input deck and the RTT Format mesh file
 *        specified therein. The addresses (i.e., opaque pointers) of both 
 *        the new CAR_CU_Interface and RTT_Format class objects are set.
 * \param self Opaque pointer to the new CAR_CU_Interface class object 
 *             (returned).
 * \param file User-input file name (must contain the RTT_format mesh file 
 *             name).
 * \param verbosity Switch used to turn detailed run-time reporting on/off.
 * \param rttFormat Opaque pointer to the new RTT_Format class object
 *                  (returned).
 */
    // Construct a CAR_CU_Interface class object from a Fortran 90 program 
    // call. This also constructs a RTT_Format class object and parses both 
    // the user-input deck and the RTT Format mesh file specified therein. 
    // The addresses (i.e., opaque pointers) of both the new CAR_CU_Interface
    // and RTT_Format class objects are set.
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

/*!
 * \brief Shadow object that destroys the CAR_CU_Interface class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 */
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
/*!
 * \brief Shadow object that returns the number of sets of grouped surface
 *        source cells for the CAR_CU_Interface class object that is 
 *        referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param ss_pos_size The number of sets of grouped surface source cells
 *                    (returned).
 */
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

/*!
 * \brief Shadow object that returns the number of grouped surface source
 *        cells in the specified set for the CAR_CU_Interface class object 
 *        that is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param surface The surface source cell set number.
 * \param ss_cells_size The number of grouped surface source cells in the 
 *                      specified set (returned).
 */
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

/*!
 * \brief Shadow object that returns the position (e.g., lox, loy, loz, hix,
 *        hiy, hiz) of all of the grouped surface source cell sets for the
 *        CAR_CU_Interface class object that is referenced by the specified 
 *        opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param pos The surface source cell set positions (returned).
 * \param ss_pos_size The string length required to return all of the surface
 *        source cell set positions.
 */
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

/*!
 * \brief Shadow object that returns the position (e.g., lox, loy, loz, hix,
 *        hiy, and hiz) of the specified set of grouped surface source cells
 *        for the CAR_CU_Interface class object that is referenced by the 
 *        specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param surface The surface source cell set number.
 * \param pos The surface source cell set position (returned).
 */
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

/*!
 * \brief Shadow object that returns the angular distribution (e.g., cosine)
 *        of the surface source cells for the CAR_CU_Interface class object 
 *        that is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param distribution The surface source angular distribution (returned).
 */
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

/*!
 * \brief Shadow object that returns the temperature of all of the grouped
 *        surface source cell sets for the CAR_CU_Interface class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param temp The surface source cell set temperatures (returned).
 * \param ss_temp_size The number of source surface cell sets.
 */
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

/*!
 * \brief Shadow object that returns the temperature of the specified grouped
 *        surface source cell set for the CAR_CU_Interface class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param surface The surface source cell set number.
 * \param temp The surface source cell set temperature (returned).
 */
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

/*!
 * \brief Shadow object that returns all of the grouped surface source cell
 *        sets for the CAR_CU_Interface class object that is referenced by
 *        the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param ss_set The surface source cell sets (returned).
 * \param ss_size The number of source cell sets.
 */
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

/*!
 * \brief Shadow object that returns the grouped surface source cells for
 *        the specified set of the CAR_CU_Interface class object that is
 *         referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param surface The surface source cell set number.
 * \param ss_set The surface source cells (returned).
 * \param ss_size The number of surface source cells.
 */
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

/*!
 * \brief Shadow object that returns the mesh volumetric sources for all of
 *        the cells in the CAR_CU_Interface class object that is referenced
 *        by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param src The cell-based volumetric sources (returned).
 * \param src_size The number of volumetric sources in the mesh (should equal
 *                 the number of cells).
 */
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

/*!
 * \brief Shadow object that returns the volumetric source for the specified
 *        cell in the CAR_CU_Interface class object that is referenced by the
 *        specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param cell The cell number.
 * \param src The cell volumetric source (returned).
 */
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

/*!
 * \brief Shadow object that returns the mesh radiation sources for all of
 *        the cells in the CAR_CU_Interface class object that is referenced
 *        by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param src The cell-based radiation sources (returned).
 * \param src_size The number of radiation sources in the mesh (should equal
 *                 the number of cells).
 */
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

/*!
 * \brief Shadow object that returns the radiation source for the specified
 *        cell in the CAR_CU_Interface class object that is referenced by the
 *        specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param cell The cell number.
 * \param src The cell radiation source (returned).
 */
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

/*!
 * \brief Shadow object that returns the cut-off time for the cell-based
 *        radiation sources in the CAR_CU_Interface class object that is
 *        referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param time The radiation source cut-off time (returned).
 */
    // Return the cut-off time for a radiation source.
    void get_car_cu_rad_s_tend_(long & self, double & time)
    {
	// Get the address of the CAR_CU_Interface class object (self).
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	time = interface->get_rad_s_tend();
    }

/*!
 * \brief Shadow object that returns the mesh analytic opacity type (e.g.,
 *        straight) for the CAR_CU_Interface class object that is referenced
 *        by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param analy_opacity The analytic opacity type (returned).
 */
    // Return the mesh analytic opacity type
    void get_car_cu_analy_opacity_(long & self, char * analy_opacity)
    {
	// Get the address of the CAR_CU Mesh (self) class object.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	string opacity = interface->get_analytic_sp_heat();

	for (int index = 0; index < opacity.size(); index++)
	{
	    * analy_opacity = opacity[index];
	    ++analy_opacity;
	}
    }

/*!
 * \brief Shadow object that returns the mesh analytic specific heat type
 *        (e.g., straight) for the CAR_CU_Interface class object that is
 *        referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param analy_spec_heat The analytic specific heat type (returned).
 */
    // Return the mesh analytic specific heat type
    void get_car_cu_analy_sp_heat_(long & self, char * analy_spec_heat)
    {
	// Get the address of the CAR_CU Mesh (self) class object.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	string spec_heat = interface->get_analytic_sp_heat();

	for (int index = 0; index < spec_heat.size(); index++)
	{
	    * analy_spec_heat = spec_heat[index];
	    ++analy_spec_heat;
	}
    }

/*!
 * \brief Shadow object that returns the mesh implicitness factor (Fleck's
 *        alpha) for the CAR_CU_Interface class object that is referenced by
 *        the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param implicit The mesh implicitness factor (returned).
 */
    // Return the mesh implicitness factor (Fleck's alpha)
    void get_car_cu_implicit_(long & self, long & implicit)
    {
	// Get the address of the CAR_CU Mesh (self) class object.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	implicit = static_cast<long>(interface->get_implicit());

    }

/*!
 * \brief Shadow object that returns the initial time step size for the
 *        CAR_CU_Interface class object that is referenced by the specified
 *        opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param time_step The initial time step size (returned).
 */
    // Return the initial time step size
    void get_car_cu_time_step_(long & self, double & time_step)
    {
	// Get the address of the CAR_CU Mesh (self) class object.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	time_step = interface->get_delta_t();

    }

/*!
 * \brief Shadow object that returns the processor capacity (cells/processor)
 *        for parallel execution for the CAR_CU_Interface class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param capacity The processor capacity (returned).
 */
    // Return the processor capacity (cells/processor).
    void get_car_cu_capacity_(long & self, long & capacity)
    {
	// Get the address of the CAR_CU Mesh (self) class object.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	capacity = static_cast<long>(interface->get_capacity());

    }

/*!
 * \brief Shadow object that returns the number of cycles to run for the
 *        CAR_CU_Interface class object that is referenced by the specified
 *        opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param num_cycles The number of cycles to run (returned).
 */
    // Return the number of cycles to run.
    void get_car_cu_num_cycles_(long & self, long & num_cycles)
    {
	// Get the address of the CAR_CU Mesh (self) class object.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	num_cycles = static_cast<long>(interface->get_max_cycle());

    }

/*!
 * \brief Shadow object that returns the print-out frequency (number of 
 *        cycles) for the CAR_CU_Interface class object that is referenced 
 *        by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Interface class object.
 * \param print_freq The print-out frequency (returned).
 */
    // Return the print frequency.
    void get_car_cu_print_freq_(long & self, long & print_freq)
    {
	// Get the address of the CAR_CU Mesh (self) class object.
	SP<CAR_CU_Interface> interface = 
	    opaque_pointers<CAR_CU_Interface>::item(self);

	print_freq = static_cast<long>(interface->get_printf());

    }

} // end extern "C"


} // end namespace rtt_amr

#endif                          // __amr_Shadow_Interface_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_Interface.cc
//---------------------------------------------------------------------------//
