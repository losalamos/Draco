//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Wed Aug 25 16:27:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for mesh library
//---------------------------------------------------------------------------//

#ifndef __mesh_Release_hh__
#define __mesh_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of mesh; 
// this can be used to get exact version information in codes that 
// use mesh
// 
//===========================================================================//

#include <string>

namespace rtt_mesh 
{
    const std::string release();
}

#endif                          // __mesh_Release_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/Release.hh
//---------------------------------------------------------------------------//
