//----------------------------------*-C++-*----------------------------------//
// Shadow_CAR_CU_Mesh_Builder.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
//---------------------------------------------------------------------------//
// @> Shadow_CAR_CU_Mesh_Builder interface file
//---------------------------------------------------------------------------//

#ifndef __mc_Shadow_CAR_CU_Mesh_Builder_cc__
#define __mc_Shadow_CAR_CU_Mesh_Builder_cc__

#include "CAR_CU_Builder.hh"
#include "CAR_CU_Interface.hh"
#include "RTT_Format.hh"
#include "CAR_CU_Mesh.hh"
#include "Shadow_Opaque_Pointers.hh"
#include <iostream>

//===========================================================================//
// Shadow_CAR_CU_Mesh_Builder - 
//
// Purpose : Provides flat interface functions to the Continuous Adaptive 
// Refinement Cartesion Unstructured Mesh Builder Class for use with Fortran.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_mc 
{
using std::cout;

// draco components
using dsxx::SP;
using rtt_imc::CAR_CU_Interface;
using rtt_format::RTT_Format;

extern "C" 
{
//---------------------------------------------------------------------------//
// CAR_CU_Builder F90 to C++ flat interface functions
//---------------------------------------------------------------------------//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
//

    // Construct a CAR_CU_Builder class from a Fortran 90 program call.  This 
    // also constructs the Coord_sys, Layout, and CAR_CU_Mesh class objects.
    // The addresses of both the new CAR_CU_Builder and CAR_CU_Mesh class 
    // objects are set. The CAR_CU_Mesh class contains member functions to 
    // return the addresses of the Coord and Layout class objects, if needed.
    void construct_car_cu_mesh_builder_(long & self, long & itf_ptr, 
					long & rttf_ptr, long & verbosity, 
					long & mesh_ptr)
    {
	bool verbose = verbosity;
	SP<CAR_CU_Interface> interface;
	SP<RTT_Format> rttFormat;
	SP<CAR_CU_Mesh> mesh;

	// Get the addresses of the CAR_CU_Interface (int_ptr) and RTT_Format 
        // (rttf_ptr) class objects.
	interface = opaque_pointers<CAR_CU_Interface>::item(itf_ptr);
	rttFormat = opaque_pointers<RTT_Format>::item(rttf_ptr);

	// Construct a new CAR_CU_Builder class object.
	SP<CAR_CU_Builder> builder(new CAR_CU_Builder(interface));

	// Construct a new CAR_CU_Mesh class object and build the mesh.
	mesh = builder->build_Mesh(rttFormat);
	interface->set_defined_surcells(builder->get_defined_surcells());
	if (verbose)
	    cout << " ** Built mesh ** " << endl;

	// return the addresses of the new CAR_CU_Builder (self) and 
	// CAR_CU_Mesh (mesh_ptr) objects.
	self = opaque_pointers<CAR_CU_Builder>::insert(builder);
	mesh_ptr = opaque_pointers<CAR_CU_Mesh>::insert(mesh);

    }

    // Destroy a CAR_CU_Builder class object from a Fortran 90 program call.
    void destruct_car_cu_mesh_builder_(long & self)
    {
	// Get the address of the CAR_CU_Builder class object (self).
	SP<CAR_CU_Builder> builder = 
	    opaque_pointers<CAR_CU_Builder>::item(self);

	// destroy the CAR_CU_Builder class object by assigning this SP to 
        // a null SP
	builder = SP<CAR_CU_Builder>();
	Ensure (!builder);

	// remove the opaque pointer to the CAR_CU_Interface class object.
	opaque_pointers<CAR_CU_Builder>::erase(self);
    }

}  // end extern "C"


}  // end namespace rtt_mc

#endif                          // __mc_Shadow_CAR_CU_Mesh_Builder_cc__

//---------------------------------------------------------------------------//
//                        end of amr_mesh/Shadow_CAR_CU_Mesh_Builder.cc
//---------------------------------------------------------------------------//
