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
// parse input

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


} // end extern "C"


} // end namespace rtt_imc

#endif                          // __mc_Shadow_CAR_CU_Interface_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_CAR_CU_Interface.cc
//---------------------------------------------------------------------------//
