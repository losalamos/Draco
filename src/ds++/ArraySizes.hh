//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ArraySizes.hh
 * \author Giovanni Bavestrelli
 * \date   Mon Apr 21 16:00:24 MDT 2003
 * \brief  A Class Template for N-Dimensional Generic Resizable Arrays.
 *
 * \sa C/C++ Users Journal, December 2000, http://www.cuj.com.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//
#ifndef rtt_dsxx_ArraySizes_hh
#define rtt_dsxx_ArraySizes_hh

#include <algorithm>

namespace rtt_dsxx
{

//! Forward declaration needed for friend declarations
template< typename T, unsigned N >
class Array;

//===========================================================================//
/*!
 * \class ArraySize
 * \brief Class that encapsulates a const unsigned (&)[N]
 */
//===========================================================================//
template< unsigned N >
class ArraySize
{
    typedef const unsigned (&UIntArrayN)[N];
    unsigned m_Dimensions[N];
    ArraySize( unsigned const (&Dimensions)[N-1],
	       unsigned const dim );

  public:
    
    //! Generate an ArraySize of the next higher dimensionality.
    ArraySize< N+1 > operator () ( unsigned dim) ; 
    //! Return a reference to the dimension array.
    operator UIntArrayN () const { return m_Dimensions; }

    friend class ArraySizes;
    friend class ArraySize<N-1>;
};


//===========================================================================//
/*!
 * \class ArraySizes
 * 
 * \brief Starting point to build a const unsigned (&)[N] on the fly
 */
//===========================================================================//
class ArraySizes
{
    //! number of dimensions (number of entries for ArraySize
    unsigned m_Dimensions[1];   
    
  public:
    
    //! default constructor for ArraySizes.
    explicit ArraySizes( unsigned const dim ) { m_Dimensions[0]=dim; }
    
    //! append next size of next dimension.
    ArraySize<2> operator () ( unsigned const dim ) { 
        return ArraySize<2>( m_Dimensions, dim ); }
};

} // end rtt_dsxx namespace

#include "ArraySizes.t.hh"

#endif // rtt_dsxx_ArraySizes_hh

