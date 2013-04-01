//----------------------------------*-C++-*----------------------------------//
/*! \file  Sync.cc
 * \author Maurice LeBrun
 * \date   Wed Jan 25 16:04:40 1995
 * \brief  Classes for forcing a global sync at the head and/or tail of a block.
 * \note   Copyright (C) 1995-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include "Sync.hh"
#include "C4_Functions.hh"  

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// HSync constructor.  Forces a global sync before continuing.
//---------------------------------------------------------------------------//

HSync::HSync( int s /*=1*/ )
{
    if (s) global_barrier();
}

//---------------------------------------------------------------------------//
// TSync destructor.  Forces a global sync before continuing. 
//---------------------------------------------------------------------------//

TSync::~TSync()
{
    if (sync) global_barrier();
}

} // end of rtt_c4

//---------------------------------------------------------------------------//
// end of Sync.cc
//---------------------------------------------------------------------------//
