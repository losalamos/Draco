//----------------------------------*-C++-*----------------------------------//
// Sync.hh
// Maurice LeBrun
// Wed Jan 25 16:04:40 1995
//---------------------------------------------------------------------------//
// @> Classes for forcing a global sync at the head and/or tail of a block.
//---------------------------------------------------------------------------//

#ifndef __c4_Sync_hh__
#define __c4_Sync_hh__

#include "c4/global.hh"

//===========================================================================//
// class HSync - Head synchronizing

// Synchronizes processes at the head of a block by doing a global sync in
// the ctor.
//===========================================================================//

class HSync {

    HSync( const HSync& );
    HSync& operator=( const HSync& );

  public:
    HSync( int s =1 );
};

//===========================================================================//
// class TSync - Tail synchronizing

// Synchronizes processes at the tail of a block by doing a global sync in
// the dtor.
//===========================================================================//

class TSync {

    TSync( const TSync& );
    TSync& operator=( const TSync& );

    int sync;

  public:
    TSync( int s =1 ) : sync(s) {}
    ~TSync();
};

//===========================================================================//
// class HTSync - Head & tail synchronizing

// Synchronizes processes at the head and tail of a block by doing a global
// sync in the ctor/dtor.
//===========================================================================//

class HTSync: public HSync, public TSync {

    HTSync( const HTSync& );
    HTSync& operator=( const HTSync& );

  public:
    HTSync( int s =1 ) : HSync(s), TSync(s) {}
};

#endif                          // __c4_Sync_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Sync.hh
//---------------------------------------------------------------------------//
