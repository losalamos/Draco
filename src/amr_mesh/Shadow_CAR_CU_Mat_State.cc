//----------------------------------*-C++-*----------------------------------//
// Shadow_CAR_CU_Mat_State.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
//---------------------------------------------------------------------------//
// @> Shadow_CAR_CU_Mat_State interface file
//---------------------------------------------------------------------------//

#ifndef __imc_Shadow_CAR_CU_Mat_State_cc__
#define __imc_Shadow_CAR_CU_Mat_State_cc__

#include "imc/Mat_State.hh"
#include "ds++/SP.hh"
#include "Shadow_Opaque_Pointers.hh"
#include "ds++/Assert.hh"
#include "CAR_CU_Mesh.hh"
#include <iostream>

//===========================================================================//
// Shadow_CAR_CU_Mat_State - 
//
// Purpose : Provides shadow interface functions to the Implicit Monte Carlo
// (imc) Mat_State Class for use with Fortran 90 codes. Note that the class 
// constructor is not shadowed herein because the Mat_State class object is 
// constructed by the Opacity_Builder class object (which is shadowed).
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_imc 
{
// draco components
using dsxx::SP;
using rtt_amr::CAR_CU_Mesh;
using std::string;

extern "C" 
{
//===========================================================================//
// Mat_State F90 to C++ flat interface functions
//===========================================================================//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
    // Destroy an Mat_State class object from a Fortran 90 program call.
    void destruct_car_cu_mat_state_(long & self)
    {
	// Get the address of the Mat_State class object (self).
	SP<Mat_State<CAR_CU_Mesh> >  material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);

	// destroy the Mat_State class object by assigning this SP to 
        // a null SP
	material = SP<Mat_State<CAR_CU_Mesh> >();
	Ensure (!material);

	// remove the opaque pointer to the Mat_State class object.
	opaque_pointers<Mat_State<CAR_CU_Mesh> >::erase(self);
    }

//===========================================================================//
// General Mat_State scalar accessor functions
//===========================================================================//
    // Return the cell density
    void get_car_cu_rho_(long & self, long & mesh_ind, long & cell, 
			 double & rho)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Mat_State<CAR_CU_Mesh> > material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_rho_!");
	Insist(mesh->num_cells() == material->num_cells(),
	       "Invalid mesh passed to get_car_cu_rho_!");

	rho = material->get_rho(icell);
    }

    // Set the cell density
    void set_car_cu_rho_(long & self, long & mesh_ind, long & cell, 
			 double & rho)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Mat_State<CAR_CU_Mesh> > material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to set_car_cu_rho_!");
	Insist(mesh->num_cells() == material->num_cells(),
	       "Invalid mesh passed to set_car_cu_rho_!");

	material->get_rho(icell) = rho;
    }

    // Return the cell temperature
    void get_car_cu_temp_(long & self, long & mesh_ind, long & cell, 
			  double & temp)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Mat_State<CAR_CU_Mesh> > material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_temp_!");
	Insist(mesh->num_cells() == material->num_cells(),
	       "Invalid mesh passed to get_car_cu_temp_!");

	temp = material->get_T(icell);
    }

    // Set the cell temperature
    void set_car_cu_temp_(long & self, long & mesh_ind, long & cell, 
			  double & temp)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Mat_State<CAR_CU_Mesh> > material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to set_car_cu_temp_!");
	Insist(mesh->num_cells() == material->num_cells(),
	       "Invalid mesh passed to set_car_cu_temp_!");

	material->get_T(icell) = temp;
    }

    // Return the gradient of the cell energy with respect to time
    void get_car_cu_dedt_(long & self, long & mesh_ind, long & cell, 
			  double & dedt)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Mat_State<CAR_CU_Mesh> > material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_dedt_!");
	Insist(mesh->num_cells() == material->num_cells(),
	       "Invalid mesh passed to get_car_cu_dedt_!");

	dedt = material->get_dedt(icell);
    }

    // Set the gradient of the cell energy with respect to time
    void set_car_cu_dedt_(long & self, long & mesh_ind, long & cell, 
			  double & dedt)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Mat_State<CAR_CU_Mesh> > material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to set_car_cu_dedt_!");
	Insist(mesh->num_cells() == material->num_cells(),
	       "Invalid mesh passed to set_car_cu_dedt_!");

	material->get_dedt(icell) = dedt;
    }

    // Return the cell specific heat
    void get_car_cu_spec_heat_(long & self, long & mesh_ind, long & cell, 
			       double & spec_heat)
    {
	// Get the addresses of the Mat_State (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Mat_State<CAR_CU_Mesh> > material = 
	    opaque_pointers<Mat_State<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_spec_heat_!");
	Insist(mesh->num_cells() == material->num_cells(),
	       "Invalid mesh passed to get_car_cu_spec_heat_!");

	spec_heat = material->get_spec_heat(icell);
    }

} // end extern "C"

} // end namespace rtt_imc

#endif                          // __imc_Shadow_CAR_CU_Mat_State_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_CAR_CU_Mat_State.cc
//---------------------------------------------------------------------------//
