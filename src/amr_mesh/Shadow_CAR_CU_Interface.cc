//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Mesh_Shadow.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
//---------------------------------------------------------------------------//
// @> CAR_CU_Mesh_Shadow interface file
//---------------------------------------------------------------------------//

#ifndef __mc_CAR_CU_Shadow_C_cc__
#define __mc_CAR_CU_Shadow_C_cc__

#include "CAR_CU_Interface.hh"
#include "RTT_Format.hh"
#include "CAR_CU_Builder.hh"
#include "CAR_CU_Mesh.hh"
#include <iostream>

//===========================================================================//
// CAR_CU_Shadow_C - 
//
// Purpose : Provides shadow interface functions to the Continuous Adaptive 
// Refinement Cartesion Unstructured Mesh Class for use with Fortran 90 codes.
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
using dsxx::Destroy;
using rtt_imc::CAR_CU_Interface;
using rtt_format::RTT_Format;

extern "C" 
{
//---------------------------------------------------------------------------//
// CAR_CU_Interface F90 to C++ shadow interface functions
//---------------------------------------------------------------------------//
// parse input for Source_Init

    // Construct a CAR_CU_Interace class from a Fortran 90 program call. This
    // also constructs a RTT_Format class object and parses both the input
    // deck and the RTT Format mesh file specified therein. The addresses of
    // both the new CAR_CU_Interace and RTT_Format class objects are set.
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

	// set the addresses of the self-referring smart pointers to the new 
	// CAR_CU_Interface (self) and RTT_Format class (rttFormat) objects.
	self = reinterpret_cast<long>(& (* interface));
	rttFormat = reinterpret_cast<long>(& (* rttMesh));

    }

    // Destruct a CAR_CU_Interace class from a Fortran 90 program call.
    void destruct_car_cu_interface_(long & self)
    {
	// Get the address of the SP to the CAR_CU_Interface class object.
	CAR_CU_Interface * interface = 
	    reinterpret_cast<CAR_CU_Interface * >(self);

	// destruct the CAR_CU_Interface class object by assigning this SP 
	// to a null SP
	delete interface;
    }


}


} // end namespace rtt_mc

#endif                          // __mc_CAR_CU_Shadow_C_cc__

//---------------------------------------------------------------------------//
//                              end of mc/CAR_CU_Shadow_C.cc
//---------------------------------------------------------------------------//
