//----------------------------------*-C++-*----------------------------------//
// Shadow_Opacity.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
/*! 
 * \file   amr_mesh/Shadow_Opacity.cc
 * \author B.T. Adams
 * \date   Mon 27 Sep 10:33:26 1999
 * \brief  Provides the C++ side of the shadow object interface functions to
 *         the Implicit Monte Carlo (imc) Opacity class for use with Fortran
 *         90 codes. The complimentary Fortran 90 shadow object interface 
 *         functions that reference the functions herein are provided in 
 *         Shadow_Opacity.f90. Note that the class constructor is not shadowed
 *         herein because the Opacity class object is constructed by the
 *         Opacity_Builder class object (which is also shadowed). An example
 *         code is also provide via link to illustrate the usage of all of the
 *         shadow object interface functions to the amr_mesh package from a
 *         Fortran 90 code.
 *
 *\sa Mark G. Gray, Randy M. Roberts, and Tom Evans, Scientific Programming,
 *   "Shadow-Object Interface Between Fortran 95 and C++", March-April 1999.
 */
//---------------------------------------------------------------------------//
// @> Shadow_Opacity interface file
//---------------------------------------------------------------------------//

#ifndef __imc_Shadow_Opacity_cc__
#define __imc_Shadow_Opacity_cc__

#include "imc/Opacity.hh"
#include "ds++/SP.hh"
#include "ds++/opaquePointers.hh"
#include "ds++/Assert.hh"
#include "Mesh.hh"
#include <iostream>

//===========================================================================//
// Shadow_Opacity - 
//
// Purpose : Provides shadow interface functions to the Implicit Monte Carlo
// (imc) Opacity Class for use with Fortran 90 codes. Note that the class 
// constructor is not shadowed because the Opacity class object is constructed
// by the Opacity_Builder class object (which is shadowed).
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_imc 
{
// draco components
using rtt_dsxx::SP;
using rtt_amr::CAR_CU_Mesh;
using rtt_dsxx::opaque_pointers;

extern "C" 
{
//===========================================================================//
// Opacity F90 to C++ flat interface functions
//===========================================================================//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
/*!
 * \brief Shadow object that destroys the Opacity class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the Opacity class object.
 */
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
/*!
 * \brief Shadow object that returns the macroscopic kappa cross section
 *        in the specified cell in the Opacity class object that is referenced
 *        by the specified opaque pointer. 
 * \param self Opaque pointer to the Opacity class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Opacity class object.
 * \param cell The cell number.
 * \param sigma_abs The cell macroscopic kappa cross section (returned).
 */
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

/*!
 * \brief Shadow object that returns the macroscopic kappa-thomson cross
 *        section in the specified cell in the Opacity class object that 
 *        is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the Opacity class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Opacity class object.
 * \param cell The cell number.
 * \param sigma_thomson The cell macroscopic kappa-thomson cross section
 *                      (returned).
 */
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

/*!
 * \brief Shadow object that returns the Planckian opacity in the specified
 *        cell in the Opacity class object that is referenced by the specified
 *        opaque pointer. 
 * \param self Opaque pointer to the Opacity class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Opacity class object.
 * \param cell The cell number.
 * \param planck The cell Planckian opacity (returned).
 */
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

/*!
 * \brief Shadow object that returns the Fleck factor in the specified
 *        cell in the Opacity class object that is referenced by the specified
 *        opaque pointer. 
 * \param self Opaque pointer to the Opacity class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Opacity class object.
 * \param cell The cell number.
 * \param fleck The cell Fleck factor (returned).
 */
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
/*!
 * \brief Shadow object that returns the Fleck factor multiplied by the 
 *        Planckian opacity in the specified cell in the Opacity class object
 *        that is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the Opacity class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Opacity class object.
 * \param cell The cell number.
 * \param fleck_planck The cell Fleck factor * Planckian opacity (returned).
 */
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

/*!
 * \brief Shadow object that returns the Fleck effective macroscopic scattering
 *        cross section in the specified cell in the Opacity class object
 *        that is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the Opacity class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Opacity class object.
 * \param cell The cell number.
 * \param sigeffscat The cell Fleck effective macroscopic scattering cross 
 *                   section (returned).
 */
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

/*!
 * \brief Shadow object that returns the Fleck effective macroscopic 
 *        absorption cross section in the specified cell in the Opacity 
 *        class object that is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the Opacity class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Opacity class object.
 * \param cell The cell number.
 * \param sigeffabs The cell Fleck effective macroscopic absorption cross 
 *                  section (returned).
 */
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

#endif                          // __imc_Shadow_Opacity_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_Opacity.cc
//---------------------------------------------------------------------------//
