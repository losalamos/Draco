//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Thomas M. Evans
// Wed Apr 14 15:36:52 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for imc library
//---------------------------------------------------------------------------//

#ifndef __imc_Release_hh__
#define __imc_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of imc; 
// this can be used to get exact version information in codes that 
// use imc
// 
//===========================================================================//

#include <string>

namespace rtt_imc 
{
    const std::string release();
}

#endif                          // __imc_Release_hh__

//---------------------------------------------------------------------------//
//                              end of Release.hh
//---------------------------------------------------------------------------//
