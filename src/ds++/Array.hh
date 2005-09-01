//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   draco/src/ds++/Array.hh
 * \author Giovanni Bavestrelli
 * \date   Mon Apr 21 16:00:24 MDT 2003
 * \brief  A Class Template for N-Dimensional Generic Resizable Arrays.
 *
 * \sa C/C++ Users Journal, December 2000, http://www.cuj.com.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//
#ifndef rtt_dsxx_Array_hh
#define rtt_dsxx_Array_hh

#include <algorithm>
#include "RefArray.hh"
#include "Assert.hh"

namespace rtt_dsxx
{

//==============================================================================
/*!
 * \class Array
 * \brief A class template for N-dimensional generic resizable arrays.
 *
 * Class Template for a Generic Resizable N Dimensional Array (for N>=2)
 * By Giovanni Bavestrelli                  
 * Copyright 1999 Giovanni Bavestrelli 
 * Any feedback is welcome, you can contact me at gbavestrelli@yahoo.com
 *
 * This is the full implementation of my array, up to the new C++ standard.
 * It uses partial specialization, and it will not work with Visual C++.
 * If you use Visual C++, take the other version of the classes.
 *
 * Example usage:
 * \code
 *     // Array sizes
 *     unisgned int Sizes3[]={1,2,3];
 *
 *     // Define some arrays.
 *     Array<int,2> const A2;
 *     Array<MyType,3> A3(Sizes3);
 *     Array<double,5> A5(ArraySizes(5)(4)(3)(2)(1));
 *
 *     // Traverse the array a'la STL and fill it.
 *     int k(0);
 *     for( Array<double,5>::iterator it=A5.begin(); it!=A5.end(); ++it)
 *     { 
 *        ++k; 
 *        *it=k/3.0;
 *     }
 *
 *     Array<double,5> CopyofA5(A5);
 *     CopyofA5.resize(ArraySizes(1)(2)(3)(4)(5)); // Resize but loose element values.
 *     CopyofA5[1][2][3][4][5] = A5[5][4][3][2][1];
 *     size_t len=A5.size();
 *
 *     // Can use STL algorithms
 *     int * pMaximum = std::max_element(A5.begin(), A5.end());
 * \endcode
 *
 * \sa C/C++ Users Journal, December 2000, <a
 * href="http://www.cuj.com/documents/s=8032/cuj0012bavestre/">A Class
 * Template for N-Dimensional Generic Resizable Arrays</a>.
 */
/*!
 * \example ds++/test/tstArray.cc
 * Example usage of the Array class.
 */
//==============================================================================

template <typename T, unsigned N>
class Array
{
  public:

    //! Type of an element of the Array.
    typedef T         value_type;
    //! Type of a reference to an element of the Array.
    typedef T       & reference;
    //! Type of a const reference to an element of the Array
    typedef T const & const_reference;
    //! Type of a pointer to an element of the Array.
    typedef T       * pointer;
    //! Type of a pointer to a const element of the Array.
    typedef T const * const_pointer;
    //! Type of an iterator into the Array.
    typedef T       * iterator;
    //! Type of a const iterator into the Array.
    typedef T const * const_iterator;
    //! Type of the number of elements in the Array.
    typedef size_t    size_type;
    //! Type of the number of elements between two iterators into the same
    //! Array. 
    typedef ptrdiff_t difference_type;
    
    //! Give access to number of dimensions
   enum  { array_dims = N };

 private:

    // The data

    T *     m_pArrayElements; //!< Pointer to actual array elements
    size_t  m_nArrayElements; //!< Total number of array elements
    
    size_t  m_NDimensions[N];  //!< Size of the N array dimensions
    size_t  m_SubArrayLen[N];  //!< Size of each subarray

 public:

    //! Default constructor
    Array<T,N>();
   
    //! Constructor used if dimensions and init value are known.
    explicit Array<T,N>( unsigned const   (&Dimensions)[N], 
			 T        const & Init=T() );

    //! Copy constructor
    Array<T,N>( Array<T,N> const & A );
    
    //! Destructor
    ~Array<T,N>() { delete [] m_pArrayElements; }
    
    //! Indexing Array
    RefArray<T,N-1>       operator []( size_t Index );
    
    //! Indexing Constant Array
    RefArray<T,N-1> const operator []( size_t Index ) const;
    
    //! Return RefArray referencing entire Array 
    RefArray<T,N> GetRefArray(); 

    //! Return constant RefArray referencing entire Array 
    RefArray<T,N> const GetRefArray() const;
    
    //! Set the size of each array dimension
    bool resize( unsigned const   (&Dimensions)[N],
		 T        const & Init=T(), 
		 bool             PreserveElems=false);

    //! Delete the complete Array
    void clear();

    //! Assignment operator
    Array<T,N> & operator = ( Array<T,N> const & A );

    //! Returns an iterator to the first element of the Array.
    iterator       begin()       { return m_pArrayElements; }
    //! Returns a  const iterator to the first element of a const Array.
    const_iterator begin() const { return m_pArrayElements; }
    //! Returns an iterator positioned just past the last element of an Array.
    iterator       end()         { return m_pArrayElements+m_nArrayElements; }
    //! Returns a const iterator positioned just past the last element of a
    //! const Array.
    const_iterator end()   const { return m_pArrayElements+m_nArrayElements; }

    //! Some more STL-like size members
    size_t size()                 const { return m_nArrayElements; }
    
    //! Return the size of each dimension, 1 to N 
    size_t size(size_t Dim) const;
    
    //! Say if the array is empty
    bool empty()                  const { return m_nArrayElements==0; }
    
    //! Return number of dimensions
    size_t dimensions()     const { return N; } 
    
    //! Swap this array with another, a'la STL 
    void swap(Array<T,N> & A);
    
  protected:
    
    // The following are protected mainly because they are not exception-safe
    // but the way they are used in the rest of the class is exception-safe
    
    //! \brief Copy the elements of another array on this one where possible
    //!  Where not possible, initialize them to a specified value Init
    void copy( Array<T,N> const & A, 
	       T          const & Init=T() );
    
    //! Initialize all the array elements
   void initialize(const T & Init=T()) { std::fill(begin(),end(),Init); }

    //! Prefer non-member operator ==, but it needs to be a friend 
   template <typename TT, unsigned NN> 
   friend bool operator == (const Array<TT,NN> & A, const Array<TT,NN> & B);
};

//---------------------------------------------------------------------------//
// External equality operators
//---------------------------------------------------------------------------//

//! Test for equality between two arrays
template <typename T, unsigned N>
bool operator == (const Array<T,N> & A, const Array<T,N> & B)
{ return std::equal(A.m_NDimensions,A.m_NDimensions+N,B.m_NDimensions)
      && std::equal(A.begin(),A.end(),B.begin()); }

//! Test for inequality between two arrays
template <typename T, unsigned N>
bool operator != (const Array<T,N> & A, const Array<T,N> & B) 
{ return !(A==B); }

//! Not implemented, meaningless to have 0 dimensions
template <typename T> class Array<T,0> { /* empty */ };
//! Not implemented, use std::vector for one dimensional arrays
template <typename T> class Array<T,1> { /* empty */ };

} // end rtt_dsxx namespace

#include "Array.t.hh"

#endif // rtt_dsxx_Array_hh

