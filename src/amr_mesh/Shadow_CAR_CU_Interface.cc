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
using rtt_imc::CAR_CU_Interface;
using rtt_format::RTT_Format;

extern "C" 
{

    void construct_car_cu_interface_(long & self, long & infile, 
				     long & length, long & verbosity, 
				     long & rttFormat)
    {
        char * file_ptr = reinterpret_cast<char *>(infile);
	int name_length = * reinterpret_cast<int *>(length);
	bool verbose = * reinterpret_cast<int *>(verbosity);
	string file;

	for (int c = 0; c < name_length; c++)
	    file.insert(c, file_ptr + c);
	cout << file << endl;
	cout << verbose << endl;
	SP<CAR_CU_Interface> interface(new CAR_CU_Interface(file,verbose));
	SP<RTT_Format> rttMesh = interface->parser();
	if (verbose)
	    cout << " ** Read input file " << file << endl;

	self = reinterpret_cast<long>( & interface);
	rttFormat = reinterpret_cast<long>( & rttMesh);

    }


}


} // end namespace rtt_mc

#endif                          // __mc_CAR_CU_Shadow_C_cc__

//---------------------------------------------------------------------------//
//                              end of mc/CAR_CU_Shadow_C.cc
//---------------------------------------------------------------------------//
