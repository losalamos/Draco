//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Wed Aug 25 16:27:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for mesh/test programs.
//---------------------------------------------------------------------------//

#ifndef __mesh_test_Release_hh__
#define __mesh_test_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of mesh/test; 
// this can be used to get exact version information in codes that 
// use mesh
// 
//===========================================================================//

#include <string>

namespace rtt_mesh_test 
{
    const std::string release();
}

#endif                          // __mesh_test_Release_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/test/Release.hh
//---------------------------------------------------------------------------//
