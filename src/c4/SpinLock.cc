//----------------------------------*-C++-*----------------------------------//
// SpinLock.cc
// Geoffrey Furnish
// Fri Dec 16 13:29:02 1994
//---------------------------------------------------------------------------//
// @> A spin lock class.  Serializes execution of a blcok.
//---------------------------------------------------------------------------//

#include "c4/SpinLock.hh"

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Constructor.  Waits for the preceeding processor to finish before
// continuing. 
//---------------------------------------------------------------------------//

SpinLock::SpinLock( int _lock /*=1*/ )
    : lock(_lock)
{
    if (lock && node)
	C4_Recv( &trash, 0, node-1, SL_Next, 0 );
}

//---------------------------------------------------------------------------//
// Here we notify the next processor in the chain that he can proceed to
// execute the block, and we go ahead about our business.
//---------------------------------------------------------------------------//

SpinLock::~SpinLock()
{
    if (lock && node < lastnode)
	C4_Send( &trash, 0, node+1, SL_Next, 0 );
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of SpinLock.cc
//---------------------------------------------------------------------------//
