//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   draco/src/ds++/RefArray.hh
 * \author Giovanni Bavestrelli
 * \brief  A Class Template for N-Dimensional Generic Resizable Arrays.
 * \sa     C/C++ Users Journal, December 2000, http://www.cuj.com.
 * \version $Id$
 */
//---------------------------------------------------------------------------//
#ifndef rtt_dsxx_RefArray_hh
#define rtt_dsxx_RefArray_hh

#include <algorithm>
#include <cstddef>

namespace rtt_dsxx
{

//! Forward declaration needed for friend declarations
template <typename T, unsigned N> class Array;

//===========================================================================//
/*!
 * \class RefArray
 * 
 * \brief Class Template for N Dimensional SubArrays within an Array.
 *
 * A RefArray is a slice out of a full Array object.  Normally a user does
 * not need to create RefArray object but they will be returned by the []
 * operator for Array objects.  
 *
 * Speed optimization:<br>
 * In some cases a RefArray object can be created to optimize access to Array
 * data.  For example, the following code snipet uses temporary RefArray
 * object to speed up the triple loop access.
 *
 * \code
 * using rtt_dsxx::Array;
 * using rtt_dsxx::RefArray;
 * using rtt_dsxx::ArraySizes;
 *
 * Array<int,3> A3( ArraySizes(10)(20)(30) );
 *
 * for( int ix=0, ival=0; ix<A3.size(1); ++ix )
 * {
 *    RefArray<int,2> A2=A3[ix];
 *    for( int iy=0; iy<A3.size(2); ++iy )
 *    {
 *       RefArray<int,1> A1=A2[iy];
 *       for( int iz=0; iz<A3.size(3); ++iz )
 *       {
 *          A1[iz] = ++ival;
 *       }
 *    }
 * }
 * \endcode
 *
 * However, if iterator access is used instead of the bracket operator the
 * above RefArray temporaries do not improve performance:
 *
 * \code
 * int ival(0);
 * for( Array<int,3>::iterator it=A3.begin(); it!=A3.end(); ++it )
 *    *it = ++ival;
 * \endcode
 *
 * \example test/tstArray.cc  
 * Unit test for Array, RefArray and ArraySizes.
 */
// revision history:
// -----------------
// 0) original
// 
// ===========================================================================

template <typename T, unsigned N>
class RefArray
{
 public:

    //! Type of an element of a RefArray.
    typedef       T   value_type;
    //! Type of a reference to an element of a RefArray.
    typedef       T & reference;
    //! Type of a reference to an element of a const RefArray.
    typedef const T & const_reference;
    //! Type of a pointer to an element of a RefArray.
    typedef       T * pointer;
    //! Type of a pointer to an element of a const RefArray.
    typedef const T * const_pointer;
    //! Type of an iterator into a RefArray.
    typedef       T * iterator;
    //! Type of an iterator into a const RefArray.
    typedef const T * const_iterator;
    //! Type of the size of a RefArray.
    typedef size_t    size_type;
    //! Type of the number of elements between two iterators into a RefArray.
    typedef ptrdiff_t difference_type;
    
    //! Give access to number of dimensions
    enum  { array_dims = N };
    
  private:
    
    size_t const * const m_pNDimensions; //!< Array dimensions
    size_t const * const m_pSubArrayLen; //!< SubArray dimensions
    T            * const m_pElements;    //!< Point to SubArray with elements within Array
    
    //! Constructor for RefArray
    RefArray<T,N>( T * pElements, 
		   const size_t * pNDimensions, 
		   const size_t * pSubArrayLen );
       
    // Disable assignment operator
    RefArray<T,N> & operator=( RefArray<T,N> const & rhs );
    
  public:
    
    //! Return a SubArray of a SubArray.
    RefArray<T,N-1>       operator [](size_t Index);
    //! Return a const SubArray of a SubArray.
    RefArray<T,N-1> const operator [](size_t Index) const;
    
    //! Return STL-like iterator
    iterator       begin()       { return m_pElements; }
    //! Return STL-like const_iterator
    const_iterator begin() const { return m_pElements; }
    //! Return STL-like end iterator
    iterator       end()         { return m_pElements+size(); }
    //! Return STL-like const end iterator
    const_iterator end()   const { return m_pElements+size(); } 
    
    //! Return size of array
    size_t size() const { return m_pNDimensions[0] * m_pSubArrayLen[0]; }
    
    //! Return size of subdimensions
    size_t size( size_t Dim ) const;

    //! Return number of dimensions 
    size_t dimensions()  const { return N; }

 protected:

    // The following are protected mainly because they are not exception-safe
    // but the way they are used in the rest of the class is exception-safe
    
    //! \brief Copy the elements of another subarray on this one where possible
    //!  Where not possible, initialize them to a specified value Init
    void copy(const RefArray<T,N> & SA, const T & Init=T())
    {
        size_t below=std::min(size(1),SA.size(1));
        size_t above=size(1);
        
        // Copy the elements we can copy
        for (size_t i=0;i<below;i++) 
            (*this)[i].copy(SA[i],Init);
        
        // Reset the elements we can't copy
        for (size_t j=below;j<above;j++)
            (*this)[j].initialize(Init);
    }

   //! Reset all the elements
   void initialize(const T & Init=T())
   {
       std::fill(begin(),end(),Init);
   }

   //! Prefer non-member operator ==, but it needs to be a friend 
   template <typename TT, unsigned NN> 
   friend bool operator == (const RefArray<TT,NN> & A, const RefArray<TT,NN> & B);

   friend class Array<T,N>;
   friend class Array<T,N+1>; 
   friend class RefArray<T,N+1>; 
};


//===========================================================================//
/*!
 * \brief Partial Specialization for Monodimensional SubArray within an Array
 */
//-----------------------------------------------------------------------------

template <typename T>
class RefArray<T,1>
{
  public:
    
    //! Type of an element of a RefArray.
    typedef       T   value_type;
    //! Type of a reference to an element of a RefArray.
    typedef       T & reference;
    //! Type of a reference to an element of a const RefArray.
    typedef const T & const_reference;
    //! Type of a pointer to an element of a RefArray.
    typedef       T * pointer;
    //! Type of a pointer to an element of a const RefArray.
    typedef const T * const_pointer;
    //! Type of an iterator into a RefArray.
    typedef       T * iterator;
    //! Type of an iterator into a const RefArray.
    typedef const T * const_iterator;
    //! Type of the size of a RefArray.
    typedef size_t    size_type;
    //! Type of the number of elements between two iterators into a RefArray.
    typedef ptrdiff_t difference_type;
    
    //! Give access to number of dimensions
    enum  { array_dims = 1 };
    
  private:
    
    size_t const * const m_pNDimensions; //!< Array dimension
    T            * const m_pElements;    //!< Point to elements within Array
    
    //! Constructor for Specialized RefArray<T,1>
    RefArray<T,1>( T            * pElements, 
		   size_t const * pNDimensions, 
		   size_t const * pSubArrayLen );
       
    // Disable assignment operator
    RefArray<T,1> & operator=( RefArray<T,1> const & rhs );

  public:


    reference       operator [](size_t Index);
    const_reference operator [](size_t Index) const;


    //! Return an iterator to the first element of the RefArray.
    iterator       begin()       { return m_pElements; }
    //! Return an iterator to the first element of the const RefArray.
    const_iterator begin() const { return m_pElements; }
    //! Return an iterator pointing just past the last element of the RefArray.
    iterator       end()         { return m_pElements+size(); }
    //! Return an iterator pointing just past the last element of a const
    //! RefArray. 
    const_iterator end()   const { return m_pElements+size(); }

    //! Return size of array
    size_t    size()       const { return m_pNDimensions[0]; } 
    //! Return size of subdimensions
    size_t    size(size_t Dim) const;
    //! Return number of dimensions 
    size_t dimensions()    const { return 1; }

  protected:
    
    // The following are protected mainly because they are not exception-safe
    // but the way they are used in the rest of the class is exception-safe
    
    //! \brief Copy the elements of another subarray on this one where possible
    //! Where not possible, initialize them to a specified value Init.
    void copy( RefArray<T,1> const & SA,  T const & Init=T() );

    //! Reset all the elements
    void initialize( T const & Init=T()) { std::fill( begin(),end(),Init ); }

    //! Prefer non-member operator ==, but it needs to be a friend 
    template <typename TT, unsigned NN> 
    friend bool operator == (const RefArray<TT,NN> & A, const RefArray<TT,NN> & B);
    
    friend class Array<T,1>;
    friend class Array<T,2>;
    friend class RefArray<T,2>; 
};

//---------------------------------------------------------------------------//
// External equality operators
//---------------------------------------------------------------------------//

//! Test for equality between two subarrays
template <typename T, unsigned N>
bool operator == (const RefArray<T,N> & A, const RefArray<T,N> & B)
{ return std::equal(A.m_pNDimensions,A.m_pNDimensions+N,B.m_pNDimensions)
      && std::equal(A.begin(),A.end(),B.begin()); }

//! Test for inequality between two subarrays
template <typename T, unsigned N>
bool operator != (const RefArray<T,N> & A, const RefArray<T,N> & B) 
{ return !(A==B); }

} // end rtt_dsxx namespace

#endif // rtt_dsxx_RefArray_hh

//---------------------------------------------------------------------------//
// end of ds++/RefArray.hh
//---------------------------------------------------------------------------//
