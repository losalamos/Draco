//----------------------------------*-C++-*----------------------------------//
// SpinLock.cc
// Geoffrey Furnish
// Fri Dec 16 13:29:02 1994
//---------------------------------------------------------------------------//
// @> A spin lock class.  Serializes execution of a blcok.
//---------------------------------------------------------------------------//

// #include "SpinLock.hh"

// namespace rtt_c4
// {

// //---------------------------------------------------------------------------//
// // Constructor.  Waits for the preceeding processor to finish before
// // continuing. 
// //---------------------------------------------------------------------------//

// SpinLock::SpinLock( int _lock /*=1*/ )
//     : trash(0),
//       lock(_lock)
// {
//     if (lock && node)
// 	receive( &trash, 0, node-1, SL_Next );
// }

// //---------------------------------------------------------------------------//
// // Here we notify the next processor in the chain that he can proceed to
// // execute the block, and we go ahead about our business.
// //---------------------------------------------------------------------------//

// SpinLock::~SpinLock()
// {
//     if (lock && node < lastnode)
// 	send( &trash, 0, node+1, SL_Next );
// }

// } // end of rtt_c4

//---------------------------------------------------------------------------//
//  end of SpinLock.cc
//---------------------------------------------------------------------------//
