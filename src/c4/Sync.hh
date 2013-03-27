//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sync.hh
 * \author Maurice LeBrun
 * \date   Wed Jan 25 16:04:40 1995
 * \brief  Classes for forcing a global sync at the head and/or tail of a block.
 * \note   Copyright (C) 1995-2013 Los Alamos National Security, LLC.
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
    TSync( const TSync& );
    TSync& operator=( const TSync& );

    int sync;

  public:
    explicit TSync( int s =1 ) : sync(s) {}
    virtual ~TSync(void);
};


// kt - Defninition of HTSync commented out.
// ------------------------------------------
// I couldn't find where this class was used.  It looks like HTSpinLock simply
// includes both TSync and HSync instead of using HTSync.

//===========================================================================//
// class HTSync - Head & tail synchronizing

// Synchronizes processes at the head and tail of a block by doing a global
// sync in the ctor/dtor.
//===========================================================================//

// class HTSync: public HSync, public TSync {

//     HTSync( const HTSync& );
//     HTSync& operator=( const HTSync& );

//   public:
//     HTSync( int s =1 ) : HSync(s), TSync(s) {}
// };

} // end of rtt_c4

#endif // __c4_Sync_hh__

//---------------------------------------------------------------------------//
// end of c4/Sync.hh
//---------------------------------------------------------------------------//
