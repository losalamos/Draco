//----------------------------------*-C++-*----------------------------------//
// Shadow_Mat_State.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
/*! 
 * \file   amr_mesh/Shadow_Mat_State.cc
 * \author B.T. Adams
 * \date   Mon 27 Sep 10:33:26 1999
 * \brief  Provides the C++ side of the shadow object interface functions to
 *         the Implicit Monte Carlo (imc) Mat_State class for use with Fortran
 *         90 codes. The complimentary Fortran 90 shadow object interface 
 *         functions that reference the functions herein are provided in 
 *         Shadow_Mat_State.f90. Note that the class constructor is not 
 *         shadowed herein because the Mat_State class object is constructed 
 *         by the Opacity_Builder class object (which is shadowed). An example
 *         code is also provide via link to illustrate the usage of all of the
 *         shadow object interface functions to the amr_mesh package from a
 *         Fortran 90 code.
 *
 *\sa Mark G. Gray, Randy M. Roberts, and Tom Evans, Scientific Programming,
 *   "Shadow-Object Interface Between Fortran 95 and C++", March-April 1999.
 */
//---------------------------------------------------------------------------//
// @> Shadow_Mat_State interface file
//---------------------------------------------------------------------------//

#ifndef __imc_Shadow_Mat_State_cc__
#define __imc_Shadow_Mat_State_cc__

#include "imc/Mat_State.hh"
#include "ds++/SP.hh"
#include "ds++/opaquePointers.hh"
#include "ds++/Assert.hh"
#include "Mesh.hh"
#include <iostream>

//===========================================================================//
// Shadow_Mat_State - 
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

/*!
 * \brief RTT implicit monte carlo (imc) namespace.
 *
 * Provides namespace protection for the Draco (RTT) implicit monte carlo
 *  utilities.
 *
 *\sa The rtt_imc namespace is documented herein only to complete the adaptive
 *    mesh refinement (amr) doxygen documentation for the mesh class. Limited
 *    shadow object interfacing of some imc classes was provided within the
 *    draco/src/amr_mesh directories to add needed functionality to the mesh
 *    class.
 */
namespace rtt_imc 
{
// draco components
using rtt_dsxx::SP;
using rtt_amr::CAR_CU_Mesh;
using rtt_dsxx::opaque_pointers;
using std::string;

extern "C" 
{
//===========================================================================//
// Mat_State F90 to C++ flat interface functions
//===========================================================================//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
/*!
 * \brief Shadow object that destroys the Mat_State class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the Mat_State class object.
 */
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
/*!
 * \brief Shadow object that returns the density in the specified cell in 
 *        the Mat_State class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the Mat_State class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Mat_State class object.
 * \param cell The cell number.
 * \param rho The cell density (returned).
 */
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

/*!
 * \brief Shadow object that sets the density in the specified cell in 
 *        the Mat_State class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the Mat_State class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Mat_State class object.
 * \param cell The cell number.
 * \param rho The cell density (supplied).
 */
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

/*!
 * \brief Shadow object that returns the temperature in the specified cell in 
 *        the Mat_State class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the Mat_State class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Mat_State class object.
 * \param cell The cell number.
 * \param temp The cell temperature (returned).
 */
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

/*!
 * \brief Shadow object that sets the temperature in the specified cell in 
 *        the Mat_State class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the Mat_State class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Mat_State class object.
 * \param cell The cell number.
 * \param temp The cell temperature (supplied).
 */
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

/*!
 * \brief Shadow object that returns the derivative of energy with respect to
 *        time in the specified cell in the Mat_State class object that is 
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the Mat_State class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Mat_State class object.
 * \param cell The cell number.
 * \param dedt The derivative of energy with respect to time in the cell 
 *             (returned).
 */
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

/*!
 * \brief Shadow object that sets the derivative of energy with respect to
 *        time in the specified cell in the Mat_State class object that is
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the Mat_State class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Mat_State class object.
 * \param cell The cell number.
 * \param dedt The derivative of energy with respect to time in the cell 
 *             (supplied).
 */
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

/*!
 * \brief Shadow object that returns the specific heat in the specified cell
 *        in the Mat_State class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the Mat_State class object.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that is
 *                 associated with this Mat_State class object.
 * \param cell The cell number.
 * \param spec_heat The cell specific heat (returned).
 */
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

#endif                          // __imc_Shadow_Mat_State_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_Mat_State.cc
//---------------------------------------------------------------------------//
