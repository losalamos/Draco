//----------------------------------*-C++-*----------------------------------//
// Shadow_RTT_Format.cc
// B.T. Adams (bta@lanl.gov)
// 4 Oct 99
//---------------------------------------------------------------------------//
// @> Shadow_RTT_Format interface file
//---------------------------------------------------------------------------//

#ifndef __mc_Shadow_RTT_Format_cc__
#define __mc_Shadow_RTT_Format_cc__

#include "RTT_Format.hh"
#include "Shadow_Opaque_Pointers.hh"
#include <iostream>

//===========================================================================//
// Shadow_RTT_Format - 
//
// Purpose : Provides shadow interface functions to the Continuous Adaptive 
// Refinement Cartesion Unstructured Mesh RTT_Format class for use with 
// Fortran 90 codes. Note that only a destructor is provided because the
// CAR_CU_Interface and a CAR_CU_Build classes (which are shadowed) provide
// full functionality to build an RTT_Format class object, parse an RTT_Format
// mesh file, and build the mesh.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_format
{

// draco components
using dsxx::SP;

extern "C" 
{
//---------------------------------------------------------------------------//
// RTT_Format F90 to C++ flat interface functions
//---------------------------------------------------------------------------//

    // Destroy a RTT_Format class object from a Fortran 90 program call.
    void destruct_car_cu_rtt_format_(long & self)
    {
	// Get the address of the RTT_Format class object (self).
	SP<RTT_Format> rttMesh = 
	    opaque_pointers<RTT_Format>::item(self);

	// destroy the RTT_Format class object by assigning this SP to 
        // a null SP
	rttMesh = SP<RTT_Format>();
	Ensure (!rttMesh);

	// remove the opaque pointer to the RTT_Format class object.
	opaque_pointers<RTT_Format>::erase(self);
    }

} // end extern "C"


} // end namespace rtt_format

#endif                          // __mc_Shadow_RTT_Format_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_RTT_Format.cc
//---------------------------------------------------------------------------//
