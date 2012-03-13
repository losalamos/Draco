//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/isFinite_pt.cc
 * \brief  Explicit template instatiations for class isFinite.
 * \note   Copyright (C) 2006-2012 Los Alamos National Security, LLC
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#include "isFinite.hh"
#include <vector>

namespace rtt_dsxx
{

template DLL_PUBLIC 
bool isInfinity( float const & x );

template DLL_PUBLIC
bool isInfinity( double const & x );

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of isFinite_pt.cc
//---------------------------------------------------------------------------//
