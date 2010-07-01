//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    ds++/destroy.hh
 * \author  Randy M. Roberts
 * \date    Thu May 20 09:23:00 1999
 * \version $Id$
 * \note   Copyright (C) 2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//

#ifndef __ds_destroy_hh__
#define __ds_destroy_hh__

#include <iterator>

namespace rtt_dsxx
{
 
//===========================================================================//
/*!
 * \fn    Destroy
 * \brief template free functions Destroy that replace the versions that are
 * no longer in the STL. 
 */
//===========================================================================//

template <class T>
inline void Destroy(T* pointer)
{
    pointer->~T();
    return;
}

template <class ForwardIterator>
inline void Destroy(ForwardIterator first, ForwardIterator last) 
{ 
    for(; first != last; ++first) 
        Destroy(&*first);
    return;
}

} // end namespace rtt_dsxx

#endif // __ds_destroy_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/destroy.hh
//---------------------------------------------------------------------------//
