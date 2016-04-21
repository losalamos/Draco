//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   ds++/Release.hh
 * \author Thomas Evans
 * \date   Thu Jul 15 09:31:44 1999
 * \brief  Header file for ds++ library release function.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//===========================================================================//
// $Id$
//===========================================================================//

#ifndef rtt_ds_Release_hh
#define rtt_ds_Release_hh

#include <string>
#include "ds++/config.h"

//===========================================================================//
/*!
 * \namespace rtt_dsxx
 *
 * \brief Namespace that contains the ds++ package classes and variables.
 *
 */
//===========================================================================//

namespace rtt_dsxx 
{

//! Query package for the release number.
DLL_PUBLIC_dsxx  const std::string release();
//! Return a list of Draco authors
DLL_PUBLIC_dsxx  const std::string author_list();
//! Return a list of Draco authors
DLL_PUBLIC_dsxx  const std::string copyright();

} // end of rtt_ds++

//! This version can be called by Fortran and wraps the C++ version.
extern "C" DLL_PUBLIC_dsxx  void ec_release( char * release_string, size_t maxlen );

#endif // rtt_ds_Release_hh

//---------------------------------------------------------------------------//
// end of ds++/Release.hh
//---------------------------------------------------------------------------//
