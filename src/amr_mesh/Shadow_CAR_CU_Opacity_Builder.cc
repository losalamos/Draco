//----------------------------------*-C++-*----------------------------------//
// Shadow_CAR_CU_Opacity_Builder.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
//---------------------------------------------------------------------------//
// @> Shadow_Opacity_Builder interface file
//---------------------------------------------------------------------------//

#ifndef __imc_Shadow_CAR_CU_Opacity_Builder_cc__
#define __imc_Shadow_CAR_CU_Opacity_Builder_cc__

#include "imc/Opacity_Builder.t.hh"
#include "CAR_CU_Interface.hh"
#include "Shadow_Opaque_Pointers.hh"
#include <iostream>

//===========================================================================//
// Shadow_CAR_CU_Opacity_Builder - 
//
// Purpose : Provides flat interface functions to the Implicit Monte Carlo
// (imc) Opacity_Builder Class for use with Fortran.
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
using rtt_amr::CAR_CU_Mesh;
using rtt_amr::CAR_CU_Interface;
using rtt_shadow::opaque_pointers;

extern "C" 
{
//---------------------------------------------------------------------------//
// Opacity_Builder F90 to C++ flat interface functions
//---------------------------------------------------------------------------//
// parse input

    // Construct a Opacity_Builder class from a Fortran 90 program call.  This 
    // also constructs the Mat_State and Opacity class objects. The addresses
    // of the new Opacity_Builder, Mat_State, and Opacity class objects are 
    // set.
    void construct_car_cu_opac_builder_(long & self, long & itf_ptr, 
					long & mesh_ptr, long & verbosity,
					long & matl_ptr, long & opac_ptr)
    {
	bool verbose = verbosity;
	SP<CAR_CU_Interface> interface;
	SP<CAR_CU_Mesh> mesh;

	// Get the addresses of the CAR_CU_Interface (int_ptr) and CAR_CU_Mesh 
        // (mesh_ptr) class objects.
	interface = opaque_pointers<CAR_CU_Interface>::item(itf_ptr);
	mesh = opaque_pointers<CAR_CU_Mesh>::item(mesh_ptr);

	// Construct a new Opacity_Builder class object.
	SP<Opacity_Builder<CAR_CU_Mesh> > 
	    opacity_builder(new Opacity_Builder<CAR_CU_Mesh>(interface));

	// Build the material state
	SP<Mat_State<CAR_CU_Mesh> > matl_state = 
	    opacity_builder->build_Mat(mesh);

	// Now build the opacity
	SP<Opacity<CAR_CU_Mesh> > opacity = 
	    opacity_builder->build_Opacity(mesh, matl_state);	
	if (verbose)
	    cout << " ** Built Mat_State and Opacity" << endl;

	// return the addresses of the new Opacity_Builder (self), Mat_State
	// (matl_ptr), and Opacity (opac_ptr) class objects.
	self = opaque_pointers<Opacity_Builder<CAR_CU_Mesh> >::
	    insert(opacity_builder);
	matl_ptr = 
	opaque_pointers<Mat_State<CAR_CU_Mesh> >::insert(matl_state);
	opac_ptr =  opaque_pointers<Opacity<CAR_CU_Mesh> >::insert(opacity);

    }

    // Destroy a Opacity_Builder class object from a Fortran 90 program call.
    void destruct_car_cu_opac_builder_(long & self)
    {
	// Get the address of the Opacity_Builder class object (self).
	SP<Opacity_Builder<CAR_CU_Mesh> > opacity_builder = 
	    opaque_pointers<Opacity_Builder<CAR_CU_Mesh> >::item(self);

	// destroy the Opacity_Builder class object by assigning this SP to 
        // a null SP
	opacity_builder = SP<Opacity_Builder<CAR_CU_Mesh> >();
	Ensure (!opacity_builder);

	// remove the opaque pointer to the CAR_CU_Interface class object.
	opaque_pointers<Opacity_Builder<CAR_CU_Mesh> >::erase(self);
    }


}  // end extern "C"


}  // end namespace rtt_imc

#endif                          // __imc_Shadow_CAR_CU_Opacity_Builder_cc__

//---------------------------------------------------------------------------//
//                       end of amr_mesh/Shadow_CAR_CU_Opacity_Builder.cc
//---------------------------------------------------------------------------//
