//----------------------------------*-C++-*----------------------------------//
// destroy.hh
// Randy M. Roberts
// Thu May 20 09:23:00 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __ds_destroy_hh__
#define __ds_destroy_hh__

#include <iterator>

namespace dsxx
{
 
//===========================================================================//
// template free functions Destroy - 
//
// Purpose :  To replace the versions that are no longer in the STL.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

    template <class T>
    inline void Destroy(T* pointer) { pointer->~T(); }

    template <class ForwardIterator>
    inline void Destroy(ForwardIterator first, ForwardIterator last) { 
	for(; first != last; ++first) 
	    Destroy(&*first);
    }

} // end namespace dsxx

#endif                          // __ds_destroy_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/destroy.hh
//---------------------------------------------------------------------------//
