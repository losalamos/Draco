//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/isFinite_pt.cc
 * \author Kent Budge
 * \date   Tue Feb 19 14:28:59 2008
 * \brief  Explicit template instatiations for class isFinite.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 * 
 * 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>
#include "isFinite.i.hh"

namespace rtt_dsxx
{

template
bool isInfinity( float const & x );

template
bool isInfinity( double const & x );

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
//                 end of isFinite_pt.cc
//---------------------------------------------------------------------------//
