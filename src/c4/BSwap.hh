//----------------------------------*-C++-*----------------------------------//
// BSwap.hh
// Geoffrey Furnish
// Tue Jan 17 13:49:50 1995
//---------------------------------------------------------------------------//
// @> A buffer swapper class.
//---------------------------------------------------------------------------//

#ifndef __c4_BSwap_hh__
#define __c4_BSwap_hh__

#include "c4/NodeInfo.hh"

C4_NAMESPACE_BEG

//===========================================================================//
// class BSwap - Exchange data buffer between two nodes

// This class provides a typesafe interface to the traditional send and
// receive functions found in message passing libraries.
//
// Note:  This class currently uses synchronous communication primitives.  It
// is not yet decided how we will handle typesafe asynchronous communication.
// We could either make a new class BSwapAsync, or we could add a constructor
// argument, or we could add an argument to the transaction methods, with
// suitable defaults.  Have to think about this issue some more.
//===========================================================================//

template<class T>
class BSwap : public NodeInfo {

    BSwap( const BSwap<T>& ) {}
    BSwap<T>& operator=( const BSwap<T>& ) { return *this; }

    int other_node;

    enum { BS_Hello = 32760,
	   BS_xmit };

  public:
    BSwap( int _other_node, int sync =0 );
    ~BSwap();

    void send( const T& d );
    void recv( T& d );

    void swap( T& d );
    void swap( const T& sd, T& rd );
};

C4_NAMESPACE_END

#endif                          // __c4_BSwap_hh__

//---------------------------------------------------------------------------//
//                              end of c4/BSwap.hh
//---------------------------------------------------------------------------//
