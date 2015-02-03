//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Slice.t.hh
 * \author Kent Budge
 * \date   Thu Jul  8 08:06:53 2004
 * \brief  Definitions of nontrivial methods of template class Slice.
 * \note   Copyright (C) 2004-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef container_Slice_t_hh
#define container_Slice_t_hh

#include "Slice.hh"

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*! 
 * \param i Index 
 * \return A reference to the element whose index is \c offset+stride*i in the
 * underlying container. 
 */
template<class T>
inline typename Slice<T>::reference
Slice<T>::operator[](typename Slice<T>::size_type i)
{
    Require(i<size());
    return first[stride*i];
}

//---------------------------------------------------------------------------//
/*! 
 * \param i Index 
 * \return A reference to the element whose index is \c offset+stride*i in the
 * underlying container. 
 */
template<class T>
inline typename Slice<T>::const_reference
Slice<T>::operator[](typename Slice<T>::size_type i) const
{
    Require(i<size());
    return first[stride*i];
}

//---------------------------------------------------------------------------//
/*! 
 * \return A reference to the element whose index is \c offset in the
 * underlying container. 
 */
template<class R>
inline typename Slice<R>::const_reference Slice<R>::front() const
{
    return first[0];
}

//---------------------------------------------------------------------------//
/*! 
 * \return A reference to the element whose index is \c
 * offset+stride*(size()-1) in the underlying container. 
 */
template<class R>
inline typename Slice<R>::const_reference Slice<R>::back() const
{
    return first[stride*(size()-1)];
}

} // end namespace rtt_dsxx

#endif // container_Slice_t_hh

//---------------------------------------------------------------------------//
// end of container/Slice.t.hh
//---------------------------------------------------------------------------//
