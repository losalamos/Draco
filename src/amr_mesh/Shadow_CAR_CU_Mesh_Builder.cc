//----------------------------------*-C++-*----------------------------------//
// Shadow_CAR_CU_Mesh_Builder.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
/*! 
 * \file   amr_mesh/Shadow_CAR_CU_Mesh_Builder.cc
 * \author B.T. Adams
 * \date   Mon 27 Sep 10:33:26 1999
 * \brief  Provides the C++ side of the shadow object interface functions to
 *         the Continuous Adaptive Refinement Cartesion Unstructured Mesh 
 *         Builder class for use with Fortran 90 codes. The complimentary 
 *         Fortran 90 shadow object interface functions that reference the 
 *         functions herein are provided in Shadow_CAR_CU_Mesh_Builder.f90.
 *         Only the class constructor and destructor are shadowed, as all 
 *         other class functions are invoked automatically by the shadow 
 *         object interface of the CAR_CU_Builder class constructor. An 
 *         example code is also provide to illustrate the usage of all of
 *         the shadow object interface functions to the amr_mesh package from
 *         a Fortran 90 code.
 *
 *\sa Mark G. Gray, Randy M. Roberts, and Tom Evans, Scientific Programming,
 *   "Shadow-Object Interface Between Fortran 95 and C++", March-April 1999.
 */
//---------------------------------------------------------------------------//
// @> Shadow_CAR_CU_Mesh_Builder interface file
//---------------------------------------------------------------------------//

#ifndef __amr_Shadow_CAR_CU_Mesh_Builder_cc__
#define __amr_Shadow_CAR_CU_Mesh_Builder_cc__

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

namespace rtt_amr 
{
using std::cout;

// draco components
using dsxx::SP;
using rtt_format::RTT_Format;
using rtt_shadow::opaque_pointers;

extern "C" 
{
//---------------------------------------------------------------------------//
// CAR_CU_Builder F90 to C++ flat interface functions
//---------------------------------------------------------------------------//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
//
/*!
 * \brief Shadow object that constructs a CAR_CU_Builder class object from
 *        a Fortran 90 program call. This also constructs the Coord_sys, 
 *        Layout, and CAR_CU_Mesh class objects. The addresses (i.e., opaque
 *        pointers) of both the new CAR_CU_Builder and CAR_CU_Mesh class
 *        objects are set. The CAR_CU_Mesh class contains member functions 
 *        to return the addresses of the Coord and Layout class objects, if
 *        needed.
 * \param self Opaque pointer to the new CAR_CU_Builder class object 
 *             (returned).
 * \param itf_ptr Opaque pointer to an existing, initialized CAR_CU_Interface
 *                class object.
 * \param rttf_ptr Opaque pointer to an existing, initialized RTT_Format class
 *                 object.
 * \param verbosity Switch used to turn detailed run-time reporting on/off.
 * \param mesh_ptr Opaque pointer to the new CAR_CU_Mesh class object
 *                 (returned).
 */
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
	if (verbose)
	    cout << " ** Built mesh ** " << endl;

	// return the addresses of the new CAR_CU_Builder (self) and 
	// CAR_CU_Mesh (mesh_ptr) objects.
	self = opaque_pointers<CAR_CU_Builder>::insert(builder);
	mesh_ptr = opaque_pointers<CAR_CU_Mesh>::insert(mesh);

    }

/*!
 * \brief Shadow object that destroys the CAR_CU_Builder class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Builder class object.
 */
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


}  // end namespace rtt_amr

#endif                          // __amr_Shadow_CAR_CU_Mesh_Builder_cc__

//---------------------------------------------------------------------------//
//                        end of amr_mesh/Shadow_CAR_CU_Mesh_Builder.cc
//---------------------------------------------------------------------------//
