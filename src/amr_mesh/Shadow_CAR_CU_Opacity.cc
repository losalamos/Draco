//----------------------------------*-C++-*----------------------------------//
// Shadow_CAR_CU_Opacity.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
//---------------------------------------------------------------------------//
// @> Shadow_Opacity interface file
//---------------------------------------------------------------------------//

#ifndef __imc_CAR_CU_Shadow_Opacity_cc__
#define __imc_CAR_CU_Shadow_Opacity_cc__

#include "imc/Opacity.hh"
#include "ds++/SP.hh"
#include "Shadow_Opaque_Pointers.hh"
#include "ds++/Assert.hh"
#include "CAR_CU_Mesh.hh"
#include <iostream>

//===========================================================================//
// Shadow_CAR_CU_Opacity - 
//
// Purpose : Provides shadow interface functions to the Continuous Adaptive 
// Refinement Cartesion Unstructured Opacity Class for use with Fortran 90 
// codes. Note that the class constructor is not shadowed because the Opacity
// class object is constructed by the Opacity_Builder class object.
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
using rtt_mc::CAR_CU_Mesh;

extern "C" 
{
//===========================================================================//
// Opacity F90 to C++ flat interface functions
//===========================================================================//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
    // Destroy an Opacity class object from a Fortran 90 program call.
    void destruct_car_cu_opacity_(long & self)
    {
	// Get the address of the Opacity class object (self).
	SP<Opacity<CAR_CU_Mesh> >  opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);

	// destroy the Opacity class object by assigning this SP to 
        // a null SP
	opacity = SP<Opacity<CAR_CU_Mesh> >();
	Ensure (!opacity);

	// remove the opaque pointer to the Opacity class object.
	opaque_pointers<Opacity<CAR_CU_Mesh> >::erase(self);
    }

//===========================================================================//
// General Opacity scalar accessor functions
//===========================================================================//
    // Return sigma = kappa * rho for a cell
    void get_car_cu_sigma_abs_(long & self, long & mesh_ind, long & cell, 
			       double & sigma_abs)
    {
	// Get the addresses of the Opacity (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_sigma_abs_!");
	Insist(mesh->num_cells() == opacity->num_cells(),
	       "Invalid mesh passed to get_car_cu_sigma_abs_!");

	sigma_abs = opacity->get_sigma_abs(icell);
    }

    // Return sigma_thomson = kappa_thomson * rho for a cell
    void get_car_cu_sigma_thomson_(long & self, long & mesh_ind, long & cell, 
				   double & sigma_thomson)
    {
	// Get the addresses of the Opacity (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_sigma_thomson_!");
	Insist(mesh->num_cells() == opacity->num_cells(),
	       "Invalid mesh passed to get_car_cu_sigma_thomson_!");

	sigma_thomson = opacity->get_sigma_thomson(icell);
    }

    // Return Planckian opacity
    void get_car_cu_planck_(long & self, long & mesh_ind, long & cell, 
			    double & planck)
    {
	// Get the addresses of the Opacity (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_planck_!");
	Insist(mesh->num_cells() == opacity->num_cells(),
	       "Invalid mesh passed to get_car_cu_planck_!");

	planck = opacity->get_planck(icell);
    }

    // Return Fleck factor
    void get_car_cu_fleck_(long & self, long & mesh_ind, long & cell, 
			   double & fleck)
    {
	// Get the addresses of the Opacity (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_fleck_!");
	Insist(mesh->num_cells() == opacity->num_cells(),
	       "Invalid mesh passed to get_car_cu_fleck_!");

	fleck = opacity->get_fleck(icell);
    }

//===========================================================================//
// General Opacity scalar operator functions
//===========================================================================//
    // Return Fleck factor * Planckian opacity
    void get_car_cu_fleck_planck_(long & self, long & mesh_ind, long & cell, 
				  double & fleck_planck)
    {
	// Get the addresses of the Opacity (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_fleck_planck_!");
	Insist(mesh->num_cells() == opacity->num_cells(),
	       "Invalid mesh passed to get_car_cu_fleck_planck_!");

	fleck_planck = opacity->fplanck(icell);
    }

    // Return Fleck effective scatter cross section
    void get_car_cu_sigeffscat_(long & self, long & mesh_ind, long & cell, 
				double & sigeffscat)
    {
	// Get the addresses of the Opacity (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_sigeffscat_!");
	Insist(mesh->num_cells() == opacity->num_cells(),
	       "Invalid mesh passed to get_car_cu_sigeffscat_!");

	sigeffscat = opacity->get_sigeffscat(icell);
    }

    // Return Fleck effective absorption cross section for a cell
    void get_car_cu_sigeffabs_(long & self, long & mesh_ind, long & cell, 
			       double & sigeffabs)
    {
	// Get the addresses of the Opacity (self) and CAR_CU Mesh (mesh_ind)
        // class objects.
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opaque_pointers<Opacity<CAR_CU_Mesh> >::item(self);
	SP<CAR_CU_Mesh> mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);

	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_car_cu_sigeffabs_!");
	Insist(mesh->num_cells() == opacity->num_cells(),
	       "Invalid mesh passed to get_car_cu_sigeffabs_!");

	sigeffabs = opacity->get_sigeffabs(icell);
    }

} // end extern "C"

} // end namespace rtt_imc

#endif                          // __imc_Shadow_CAR_CU_Opacity_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_CAR_CU_Opacity.cc
//---------------------------------------------------------------------------//
