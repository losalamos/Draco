//----------------------------------*-C++-*----------------------------------//
// BSwap.cc
// Geoffrey Furnish
// Tue Jan 17 13:49:51 1995
//---------------------------------------------------------------------------//
// @> A buffer swapper class.
//---------------------------------------------------------------------------//

#include "c4/BSwap.hh"
#include "c4/global.hh"

#include "ds++/Assert.hh"

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Constructor.  Set up internal state and prepare for communiction.  Sync
// with partner if requested.
//---------------------------------------------------------------------------//

template<class T>
BSwap<T>::BSwap( int _other_node, int sync /*=0*/ )
    : other_node(_other_node)
{
    if (sync) {
    // It's probably not a good idea to use this in the general case, since I
    // suppose deadlocks are probably pretty easy to trigger.  But if the
    // pattern of hookups is well understood, it may be possible to use this
    // to effect.

	int trash;

	if (node < other_node) {

	// Send salutation.
	    Send( &trash, 0, other_node, BS_Hello, group );

	// Wait for acknoledgement.
	    Recv( &trash, 0, other_node, BS_Hello, group );

	} else if ( node > other_node ) {

	// Wait for salutation.
	    Recv( &trash, 0, other_node, BS_Hello, group );

	// Send acknoledgement.
	    Send( &trash, 0, other_node, BS_Hello, group );

	} else {

	// Swap with self?  Does this guy know what he's doing?
	    throw( "Excuuuuuse me?" );
	}
    }
}

//---------------------------------------------------------------------------//
// Tear down any internal machinery.
//---------------------------------------------------------------------------//

template<class T>
BSwap<T>::~BSwap()
{
}

//---------------------------------------------------------------------------//
// Send a data buffer to the other node.
//---------------------------------------------------------------------------//

template<class T>
void BSwap<T>::send( const T& d )
{
//     C4_Send( (void *) &d, sizeof(d), other_node, BS_xmit, group );
    Send( &d, 1, other_node, BS_xmit, group );
}

//---------------------------------------------------------------------------//
// Receive a data buffer from the other node.
//---------------------------------------------------------------------------//

template<class T>
void BSwap<T>::recv( T& d )
{
// Needs attention.  Can we convert this to Recv( T *, ... ) ???
    int len = Recv( (char *) &d, sizeof(d), other_node, BS_xmit, group );

    Assert( len == sizeof(d) );
}

//---------------------------------------------------------------------------//
// Perform "in-place" exchange of data buffer with other node.
//---------------------------------------------------------------------------//

template<class T>
void BSwap<T>::swap( T& d )
{
    send( d );
    recv( d );
}

//---------------------------------------------------------------------------//
// Swap data with the other node, but don't overwrite the send buffer with the
// received data.
//---------------------------------------------------------------------------//

template<class T>
void BSwap<T>::swap( const T& sd, T& rd )
{
    send( sd );
    recv( rd );
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of BSwap.cc
//---------------------------------------------------------------------------//
