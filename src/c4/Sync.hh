//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sync.hh
 * \author Maurice LeBrun
 * \date   Wed Jan 25 16:04:40 1995
 * \brief  Classes for forcing a global sync at the head and/or tail of a block.
 * \note   Copyright (C) 1995-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_Sync_hh__
#define __c4_Sync_hh__

#include "ds++/config.h"

namespace rtt_c4
{

//===========================================================================//
// class HSync - Head synchronizing

// Synchronizes processes at the head of a block by doing a global sync in
// the ctor.
//===========================================================================//

class DLL_PUBLIC HSync
{
	// Disable copy constructor and assignment operators.
    HSync( const HSync& );
    HSync& operator=( const HSync& );

  public:
    explicit HSync( int s =1 );
    virtual ~HSync(void) {/*empty*/};
};

//===========================================================================//
// class TSync - Tail synchronizing

// Synchronizes processes at the tail of a block by doing a global sync in
// the dtor.
//===========================================================================//

class DLL_PUBLIC TSync
{
	// Disable copy constructor and assignment operators.
    TSync( const TSync& );
    TSync& operator=( const TSync& );

    int sync;

  public:
    explicit TSync( int s =1 ) : sync(s) {}
    virtual ~TSync(void);
};

} // end of rtt_c4

#endif // __c4_Sync_hh__

//---------------------------------------------------------------------------//
// end of c4/Sync.hh
//---------------------------------------------------------------------------//
