//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SpinLock.hh
 * \author Geoffrey Furnish
 * \date   Fri Dec 16 13:29:01 1994
 * \brief  A spin lock class.  Serializes execution of a blcok.
 * \note   Copyright (C) 1994-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_SpinLock_hh__
#define __c4_SpinLock_hh__

#include "NodeInfo.hh"
#include "Sync.hh"

namespace rtt_c4
{

//===========================================================================//
/*! \class SpinLock Serialize execution of a block.
 *
 * This class enables you to get a block of code to execute serially.  Each
 * processor begins to executes the block only after the one before it is
 * finished.
 */
//===========================================================================//

class DLL_PUBLIC SpinLock : public NodeInfo
{
    // Disable copy constructor and assignment operators.
    SpinLock(            const SpinLock& );
    SpinLock& operator=( const SpinLock& );

    enum { SL_Next = 92874 };

    int trash;
    int lock;

  public:
    /*! \brief Constructor.  Waits for the preceeding processor to finish
     *  before continuing. */
    SpinLock( int _lock =1 )
        : trash(0), lock(_lock) {
            if (lock && node) receive( &trash, 0, node-1, SL_Next ); }
    
    /*! \brief Destructor Here we notify the next processor in the chain that
     *   he can proceed to execute the block, and we go ahead about our
     *  business. */
    ~SpinLock() {
        if (lock && node < lastnode) send( &trash, 0, node+1, SL_Next ); }
};

//===========================================================================//
/*! \class HSyncSpinLock Serialize a block, syncing at top.
 *
 * A spinlock that forces a global sync at the head of the block.
 */
//===========================================================================//

class DLL_PUBLIC HSyncSpinLock : public HSync, public SpinLock 
{
	// disable copy and assignment operators.
    HSyncSpinLock( const HSyncSpinLock& );
    HSyncSpinLock& operator=( const HSyncSpinLock& );

  public:
    HSyncSpinLock( int l =1 ) : HSync(l), SpinLock(l) {}
};

//===========================================================================//
/*! \class TSyncSpinLock Serialize a block, syncing at bottom.
 *
 * A spinlock that forces a global sync at the tail of the block.
 */
//===========================================================================//

class DLL_PUBLIC TSyncSpinLock : public TSync, public SpinLock 
{
	// disable copy and assignment operators.
    TSyncSpinLock( const TSyncSpinLock& );
    TSyncSpinLock& operator=( const TSyncSpinLock& );

  public:
    TSyncSpinLock( int l =1 ) : TSync(l), SpinLock(l) {}
};

//===========================================================================//
/*! \class HTSyncSpinLock Serialize a block, syncing at top and bottom.
 *
 * A spinlock that forces a global sync at the head and tail of the block.
 */
//===========================================================================//

class DLL_PUBLIC HTSyncSpinLock : public HSync, public TSync, public SpinLock 
{
	// disable copy and assignment operators.
    HTSyncSpinLock( const HTSyncSpinLock& );
    HTSyncSpinLock& operator=( const HTSyncSpinLock& );

  public:
    HTSyncSpinLock( int l =1 ) : HSync(l), TSync(l), SpinLock(l) {}
};

} // end of rtt_c4

#endif // __c4_SpinLock_hh__

//---------------------------------------------------------------------------//
// end of c4/SpinLock.hh
//---------------------------------------------------------------------------//
