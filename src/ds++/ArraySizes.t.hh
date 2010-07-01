//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   draco/src/ds++/ArraySizes.t.hh
 * \author Kelly Thompson
 * \date   Mon Apr 21 16:00:24 MDT 2003
 * \brief  ArraySizes template implementation.
 * \note   Copyright (C) 2003-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_ArraySizes_t_hh
#define rtt_dsxx_ArraySizes_t_hh

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
// ArraySize implementations
//---------------------------------------------------------------------------//

template < unsigned N >
ArraySize< N >::ArraySize( unsigned const (&Dimensions)[N-1],
			   unsigned const dim )
{
    std::copy(&Dimensions[0],&Dimensions[N-1],m_Dimensions);
    m_Dimensions[N-1]=dim;
}

//---------------------------------------------------------------------------//
/*! 
 * 
 * \param dim Size associated with the new dimension.
 */

template < unsigned N >
ArraySize< N+1 > ArraySize< N >::operator () ( unsigned dim ) 
{ 
    return ArraySize< N+1 >( m_Dimensions, dim );
}

//---------------------------------------------------------------------------//
} // end namespace rtt_dsxx
//---------------------------------------------------------------------------//

#endif // rtt_dsxx_ArraySizes_t_hh

//---------------------------------------------------------------------------//
//  end of ArraySizes.t.hh
//---------------------------------------------------------------------------//
