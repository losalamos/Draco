//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   draco/src/ds++/Array.t.hh
 * \author Kelly Thompson
 * \date   Mon Apr 21 16:00:24 MDT 2003
 * \brief  Array template implementation.
 * \note   Copyright (C) 2003-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_Array_t_hh
#define rtt_dsxx_Array_t_hh

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
RefArray<T,N>::RefArray( T * pElements,
			 size_t const * pNDimensions,
			 size_t const * pSubArrayLen )
    : m_pNDimensions( pNDimensions ),
      m_pSubArrayLen( pSubArrayLen ),
      m_pElements(    pElements    )
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
RefArray<T,1>::RefArray( T            * pElements,
			 size_t const * pNDimensions,
			 size_t const * Remember(pSubArrayLen) )
    : m_pNDimensions( pNDimensions ),
      m_pElements(    pElements    )      
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
 * \param Dim An unsignedeger specifying the index of the dimension that
 * is being queried.
 * \return A size_t that specifies the length of Array for dimension Dim.
 */
template< typename T >
typename RefArray<T,1>::size_type RefArray<T,1>::size( size_t Remember(Dim) ) const
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
// Array<T,N> implementation
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Default constructor for Array
 */
template< typename T, unsigned N >
Array<T,N>::Array()
    : m_pArrayElements( NULL ),
      m_nArrayElements( 0    )
{
    std::fill( m_NDimensions,m_NDimensions+N, 0 );
    std::fill( m_SubArrayLen,m_SubArrayLen+N, 0 );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor for Array if size and init value are specified.
 * \param Dimensions An array of N values representing the size of the N
 *        dimensions.
 * \param Init Value used to initialize all members of the Array.
 */
template< typename T, unsigned N >
Array<T,N>::Array( unsigned const (&Dimensions)[N],
		   T        const & Init )
    : m_pArrayElements( NULL ),
      m_nArrayElements( 0    )
{
    std::fill( m_NDimensions,m_NDimensions+N, 0 );
    std::fill( m_SubArrayLen,m_SubArrayLen+N, 0 );
    resize( Dimensions, Init );
}

//! Copy constructor
template< typename T, unsigned N >
Array<T,N>::Array( Array<T,N> const & A )
    : m_pArrayElements( NULL ),
      m_nArrayElements( 0    )
{
    std::fill( m_NDimensions,m_NDimensions+N, 0 );
    std::fill( m_SubArrayLen,m_SubArrayLen+N, 0 );
    
    Array<T,N> Temp;

    // Note:  The following command failed on QSC:
    //    if( !A.empty() && Temp.resize(A.m_NDimensions) )
    // To fix the problem, I needed to provide an explicit promotion from
    // std::size_t (&)[N] to unsigned (&)[N].  This is done by creating
    // Temp_dimensions as an unsigned (&)[N] and copy the data from
    // A.m_NDimensions into the the new dimensions array.  Finally, the
    // temporary dimensions array is used to resize Temp.
    unsigned Temp_dimensions[N];
    std::copy(A.m_NDimensions,A.m_NDimensions+N,Temp_dimensions);
    // if( !A.empty() ) // I believe this test is redundant with the next test.
	if( Temp.resize(Temp_dimensions) )
	    std::copy(A.begin(),A.end(),Temp.begin());
    swap( Temp );
}

//! Indexing Array
template< typename T, unsigned N >
RefArray<T,N-1> Array<T,N>::operator []( size_t Index ) 
{  
    Require( m_pArrayElements );
    Require( Index < m_NDimensions[0] );
    return RefArray<T,N-1>( &m_pArrayElements[Index*m_SubArrayLen[0]],
			    m_NDimensions+1,
			    m_SubArrayLen+1 );
}

//! Indexing Constant Array
template< typename T, unsigned N >
RefArray<T,N-1> const Array<T,N>::operator []( size_t Index ) const 
{  
    Require(m_pArrayElements);
    Require(Index<m_NDimensions[0]);
    return RefArray<T,N-1>(&m_pArrayElements[Index*m_SubArrayLen[0]],
			   m_NDimensions+1,m_SubArrayLen+1);
}

//! Return RefArray referencing entire Array 
template< typename T, unsigned N >
RefArray<T,N> Array<T,N>::GetRefArray() 
{  
    Require(m_pArrayElements);
    return RefArray<T,N>(m_pArrayElements,m_NDimensions,m_SubArrayLen);
}

//! Return constant RefArray referencing entire Array 
template< typename T, unsigned N >
RefArray<T,N> const Array<T,N>::GetRefArray() const 
{  
    Require(m_pArrayElements);
    return RefArray<T,N>(m_pArrayElements,m_NDimensions,m_SubArrayLen);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Resize the Array.
 * 
 * \param Dimensions A c-style array with N entries.  Each entry spcifies the
 * length in each dimension for the new size.
 * \param Init Initialization value for new elements.
 * \param PreserveElems If true the algorithm will preserve current data.
 * \return Array with a new size.
 */
template< typename T, unsigned N >
bool Array<T,N>::resize( unsigned const   (&Dimensions)[N], 
			 T        const & Init, 
			 bool             PreserveElems )
{
    Array<T,N> Temp; 
    
    // Calculate all the information you need to use the array
    Temp.m_nArrayElements=1;
    for (size_t i=0;i<N;i++)
    {
	if (Dimensions[i]==0)
            return false; // Check that no dimension was zero 
	Temp.m_nArrayElements*=Dimensions[i]; 
	Temp.m_NDimensions[i]=Dimensions[i];
	Temp.m_SubArrayLen[i]=1;              
	for (size_t k=N-1;k>i;k--)
            Temp.m_SubArrayLen[i]*=Dimensions[k];
    }  
    
    // Allocate new elements, let exception propagate
    Temp.m_pArrayElements=new T[Temp.m_nArrayElements];
    
    // Some compilers might not throw exception if allocation fails
    // The coverage checker marks this line as not covered for
    // Temp.m_pArrayElements == false.  This is okay because we don't want
    // this check to be true for all of our compilers.
    if (!Temp.m_pArrayElements)
	return false;
    
    // Copy the elements from the previous array if requested
    // [kt 9/29/03] Previously this if statement included an additional
    // check:  "&& !empty()".  I removed this test because earlier in this
    // function we ensure that all dimensions > 0.  If we get this far, test
    // will always be true.
    if( PreserveElems ) 
	Temp.copy(*this,Init);
    // Otherwise initialize them to the specified value
    else 
	Temp.initialize(Init);
    
    // Now swap this object with the temporary
    swap(Temp);
    
    return true; 
}

//! Delete the complete Array
template< typename T, unsigned N >
void Array<T,N>::clear()
{ 
    delete [] m_pArrayElements; 
    m_pArrayElements=NULL; 
    m_nArrayElements=0;
    
    std::fill(m_NDimensions,m_NDimensions+N,0);
    std::fill(m_SubArrayLen,m_SubArrayLen+N,0);
}

//! Assignment operator
template< typename T, unsigned N >
Array<T,N> & Array<T,N>::operator = ( Array<T,N> const & A )
{
    if (&A!=this) // For efficiency  
    {
        Array<T,N> Temp(A);
        swap(Temp);
    }
    return *this;
}

//! Return the size of each dimension, 1 to N 
template< typename T, unsigned N >
typename Array<T,N>::size_type Array<T,N>::size( size_t Dim ) const 
{  
    Require(Dim>=1);
    Require(Dim<=N); 
    return m_NDimensions[Dim-1];  
}

//! Swap this array with another, a'la STL 
template< typename T, unsigned N >
void Array<T,N>::swap( Array<T,N> & A )
{
    std::swap(m_pArrayElements,A.m_pArrayElements); 
    std::swap(m_nArrayElements,A.m_nArrayElements); 

    std::swap_ranges(m_NDimensions,m_NDimensions+N,A.m_NDimensions);
    std::swap_ranges(m_SubArrayLen,m_SubArrayLen+N,A.m_SubArrayLen);
}


//---------------------------------------------------------------------------//
/*! 
 * \brief Copy the elements of another array on this one
 *
 * Copy the elements of another array on this one where possible
 *  Where not possible, initialize them to a specified value Init
 *
 * \param A Copy elements from this Array
 * \param Init Initial value for elements of current Array that cannot be
 * copied from A.
 * \return void
 */
template< typename T, unsigned N >
void Array<T,N>::copy( Array<T,N> const & A, 
		       T          const & Init )
{
    size_t below=std::min(size(1),A.size(1));
    size_t above=size(1);
    
    // Copy the elements we can copy
    for (size_t i=0;i<below;i++)
	(*this)[i].copy(A[i],Init);
    
    // Reset the elements we can't copy
    for (size_t j=below;j<above;j++)
	(*this)[j].initialize(Init);
}

//---------------------------------------------------------------------------//
} // end namespace rtt_dsxx
//---------------------------------------------------------------------------//

#endif // rtt_dsxx_Array_t_hh

//---------------------------------------------------------------------------//
//  end of Array.t.hh
//---------------------------------------------------------------------------//
