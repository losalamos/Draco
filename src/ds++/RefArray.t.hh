//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   draco/src/ds++/RefArray.t.hh
 * \author Kelly Thompson
 * \date   Mon Apr 21 16:00:24 MDT 2003
 * \brief  Array template implementation.
 * \note   Copyright (C) 2003-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_RefArray_t_hh
#define rtt_dsxx_RefArray_t_hh

#include "RefArray.hh"

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
// RefArray Implementation
//---------------------------------------------------------------------------//

/*!
 * \brief Constructor for RefArray.
 *
 * \param pElements Pointer to SubArray with elements within Array
 * \param pNDimensions Array dimensions
 * \param pSubArrayLen SubArray dimensions
 * \return RefArray
 */
template<typename T, unsigned N>
RefArray<T,N>::RefArray<T,N>( T            * pElements,
                              size_t const * pNDimensions,
                              size_t const * pSubArrayLen )
    : m_pElements(    pElements    ),
      m_pNDimensions( pNDimensions ),
      m_pSubArrayLen( pSubArrayLen )   
{
    // Ensure that no arguments are zero.
    // The RefArray must point to real data.
    Ensure( m_pElements );
    Ensure( m_pNDimensions ); 
    Ensure( m_pSubArrayLen );

    // Ensure that the RefArray is not zero length.
    Ensure( m_pNDimensions[0] > 0 );
    Ensure( m_pSubArrayLen[0] > 0 );
}

//---------------------------------------------------------------------------//

/*!
 * \brief The bracket operator will return a slice of a SubArray as a new
 * smaller SubArray.
 *
 * \param Index specify the <i>slice</i> requested from the larger SubArray.
 * \return RefArray that has a smaller dimensionality.
 */
template<typename T, unsigned N>
RefArray<T,N-1> RefArray<T,N>::operator []( size_t Index )
{
    Require( m_pElements );
    Require( Index < m_pNDimensions[0] );
    return RefArray<T,N-1>( &m_pElements[ Index*m_pSubArrayLen[ 0 ] ],
			    m_pNDimensions+1,
			    m_pSubArrayLen+1 );
}

/*!
 * \brief The const bracket operator will return a slice of a SubArray as a new
 * smaller const SubArray.
 *
 * \param Index specify the <i>slice</i> requested from the larger SubArray.
 * \return const RefArray that has a smaller dimensionality.
 */
template<typename T, unsigned N>
RefArray<T,N-1> const RefArray<T,N>::operator []( size_t Index ) const
{
    Require( m_pElements );
    Require( Index < m_pNDimensions[0] );
    return RefArray<T,N-1>( &m_pElements[ Index*m_pSubArrayLen[ 0 ] ],
			    m_pNDimensions+1,
			    m_pSubArrayLen+1 );
}

//---------------------------------------------------------------------------//

/*!
 * \brief Return the size of dimension Dim for the current SubArray.
 *
 * \param Dim Dimension of array this is to be queried.
 * \return The size of the SubArray for dimension Dim.
 */
template<typename T, unsigned N>
typename RefArray<T,N>::size_type RefArray<T,N>::size( size_t Dim ) const
{
    Require( Dim >= 1 );
    Require( Dim <= N );
    return m_pNDimensions[ Dim-1 ];
}

//---------------------------------------------------------------------------//
// Specialized RefArray<T,1>
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for RefArray.
 *
 * \param pElements Pointer to SubArray with elements within Array
 * \param pNDimensions Array dimensions
 * \param pSubArrayLen SubArray dimensions
 * \return RefArray
 */

template<typename T>
RefArray<T,1>::RefArray<T,1>( T               * pElements,
			      size_t const * pNDimensions,
			      size_t const * pSubArrayLen )
    : m_pElements(    pElements    ),
      m_pNDimensions( pNDimensions )
{
    // Ensure that no arguments are zero.
    // The RefArray must point to real data.
    Ensure( m_pElements );
    Ensure( m_pNDimensions ); 
    Ensure( pSubArrayLen );

    // Ensure that the RefArray is not zero length.
    Ensure( m_pNDimensions[0] > 0 );
    Ensure( pSubArrayLen[0] == 1 );
}

/*!
 * \brief The bracket operator will return the data.
 * smaller SubArray.
 *
 * \param Index specify the <i>element</i> requested from the larger SubArray.
 * \return reference to actual data element.
 */
template< typename T >
typename RefArray<T,1>::reference RefArray<T,1>::operator []( size_t Index )
{
    Require( m_pElements );
    Require( Index < m_pNDimensions[0] );
    return m_pElements[Index];
}

/*!
 * \brief The const bracket operator will return the data.
 *
 * \param Index specify the <i>element</i> requested from the larger SubArray.
 * \return const reference to actual data element.
 */
template< typename T >
typename RefArray<T,1>::const_reference RefArray<T,1>::operator []( size_t Index ) const
{
    Require( m_pElements );
    Require( Index < m_pNDimensions[0] );
    return m_pElements[Index];
}

/*! 
 * \brief Return the size of subdimensions.
 * 
 * \param Dim An size_t specifying the index of the dimension that
 * is being queried.
 * \return A size_type that specifies the length of Array for dimension Dim.
 */
template< typename T >
typename RefArray<T,1>::size_type RefArray<T,1>::size( size_t Dim ) const
{ 
    Require( Dim == 1 );
    return m_pNDimensions[0];
}

/*! 
 * \brief 
 * 
 * \param SA
 * \param Init 
 * \return void
 */
template< typename T >
void RefArray<T,1>::copy( RefArray<T,1> const & SA, T const & Init )
{ 
    size_t below( std::min(size(1), SA.size(1)) );
    size_t above( size(1) );
    
    // Copy the elements we can copy
    for( size_t i=0; i<below; i++ )
	m_pElements[i]=SA.m_pElements[i];
    
    // Reset the elements we can't copy
    for( size_t j=below; j<above; j++ )
	m_pElements[j]=Init;
}

//---------------------------------------------------------------------------//
} // end namespace rtt_dsxx
//---------------------------------------------------------------------------//

#endif // rtt_dsxx_RefArray_t_hh

//---------------------------------------------------------------------------//
//  end of RefArray.t.hh
//---------------------------------------------------------------------------//
