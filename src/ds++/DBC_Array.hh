//----------------------------------*-C++-*----------------------------------//
/*!
  \file    ds++/DBC_Array.hh
  \brief   Declaration of class DBC_Array
  \note    Copyright 2005-2012 Los Alamos National Security, LLC
  \version $Id$
*/
//---------------------------------------------------------------------------//

#ifndef DBC_Array_hh
#define DBC_Array_hh

#include "Assert.hh"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <limits>
#include <cstddef>

namespace rtt_dsxx
{
/*!
  Support classes and functions for DBC_Array.  This is necessary to
  implement the iterator-range constructor and assign member functions.
*/
namespace DBCA_Support
{


//---------------------------------------------------------------------------//
/*!  
  All DBC_Array constructors use this function to do an initial memory
  allocation of an array of type \c T of length \p n.  If \p init is \c true,
  then the data is copied into the newly allocated memory, otherwise the
  memory is initialized according to the rules in section 5.3.4.15 of the
  standard.

  This works differently than the allocator for a std::vector when T is
  non-POD.  Here, an array is default-initialized by array new, then a value
  is copied into each location (if init == true) using operator=().  In
  std::vector, each location is initialized using placement new and T(value).
  In a vector where no initialization value is given, the T() ctor is called
  once to create a temporary, then N calls to the T(T&) copy ctor are made,
  then one call to ~T() to kill the temporary.  This differs in that N calls
  are made to T() by array new, then N calls are made to op=(T&) by fill_n.

  \note\b Note: should you be tempted to replace new/delete[] with malloc/free,
  you had better know what you are doing with placement new and/or
  uninitialized_fill.  malloc/free works fine with PODs, but it doesn't
  construct any non-trivial classes!  When I wrote this, the performance of
  the two styles was pretty much equivalent anyway. 
*/
template<class T> void 
common_construct(T*& d_ptr, 
		 size_t& d_size,
		 const size_t n, 
		 const bool init, 
		 T const& value = T())
{
    d_size = n;
    
    if(d_size) 
    {
	d_ptr = new T[n];

	if(!d_ptr)
	{
	    std::perror("DBC_Array<T>::common_construct()");
	    Insist_ptr(d_ptr, "allocation failure");
	}
	
	if(init)
	    std::fill_n(d_ptr, n, value);
    } 
    else
    {
	d_ptr = 0;
    }
}


//---------------------------------------------------------------------------//
/*!  
  Any DBC_Array function that changes the size of the array uses this
  function to do so, creating a new array of type \c T of length \p n.  If \p
  init is \c true, then the data is copied into the newly allocated memory,
  otherwise the data is initialized according to the rules in section
  5.3.4.15 of the standard.

  This differs from the constructor in that it deletes existing memory if a
  different size is requested.  Note that there is no concept of capacity
  vs. size, as in a std::vector. 

  \sa common_construct   
*/
template<class T> void 
common_resize(T*& d_ptr, 
	      size_t& d_size,
	      const size_t n, 
	      const bool init, 
	      T const& value = T())
{
    if(d_size != n)
    {
	if(d_ptr) 
	{
	    delete[] d_ptr; 
	    d_ptr = 0;
	}
	d_size = n;
	if(d_size)
	{
	    d_ptr = new T[n];
		    
	    if(!d_ptr)
	    {
		std::perror("DBC_Array<T>::common_resize()");
		Insist_ptr(d_ptr, "allocation failure");
	    }

	} 
    }
    if(init)
	std::fill_n(d_ptr, n, value);
}


//---------------------------------------------------------------------------//
//! Implement the iterator-range constructor
template<class T, class IV, bool> struct initialize_dispatcher
{
    inline static void 
    do_it(T*& d_ptr, size_t& d_size, IV first, IV last)
    {
	common_construct<T>(d_ptr, d_size, std::distance(first, last), false);
	if(d_size)
	{
	    std::copy(first, last, d_ptr);
	}
    }
};


//---------------------------------------------------------------------------//
//! Implement the size,value constructor masquerading as iterator-range
template<class T, class IV> struct initialize_dispatcher<T,IV,true>
{
    inline static void 
    do_it(T*& d_ptr, size_t& d_size, IV n, IV val)
    {
	common_construct<T>(d_ptr, d_size, n, true, val);
    }
};


//---------------------------------------------------------------------------//
//! Implement the iterator-range assign()
template<class T, class IV, bool> struct assign_dispatcher
{
    inline static void
    do_it(T*& d_ptr, size_t& d_size, IV first, IV last)
    {
	common_resize<T>(d_ptr, d_size, std::distance(first, last), false);
	if(d_size)
	{
	    std::copy(first, last, d_ptr);
	}
    }
};


//---------------------------------------------------------------------------//
//! Implement the size,value assign() masquerading as iterator-range
template<class T, class IV> struct assign_dispatcher<T,IV,true>
{
    inline static void
    do_it(T*& d_ptr, size_t& d_size, IV n, IV val)
    {
	common_resize<T>(d_ptr, d_size, n, true, val);
    }
};

//---------------------------------------------------------------------------//
} // End of namespace rtt_dsxx::DBCA_Support 
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
//! Another "almost" replacement for std::vector
/*! 
  Why another replacement for std::vector?  I just want something simple.  No
  fancy allocators, no weird address manipulation, etc.  This is as close to
  using a malloc'd array as possible.  Some things that are different from
  std::vector: 

  - No default initialization by default (got that? ;-)
  - Different ctor behavior (see the comments above).
  - Not much dynamic length stuff (e.g. no reserve, push_back, pop_back, etc..)
  - \c clear() lets go of the memory!
  - Fast ;-)

  This class was designed to represent field data.  All constructors work
  fine with zero length allocations.
*/
template<class T> class DBC_Array
{
  public:
    typedef T                 value_type;
    typedef T &               reference;
    typedef T const &         const_reference;
    typedef T*                iterator;
    typedef T const *         const_iterator;
    typedef size_t            size_type;
    typedef ptrdiff_t         difference_type;
    typedef value_type*       pointer;
    typedef value_type const* const_pointer;
  public:
    //! Create an empty array
    DBC_Array(void)  
	: d_ptr(0), 
	  d_size(0) 
    {/* empty */}

    //! Create an uninitialized array of size \p n
    /*! 
     * If T is a POD, the contents of the array will be arbitrary.  Otherwise,
     * the default ctor is called for each array element.
     */
    explicit DBC_Array(const size_type n)
         : d_ptr(0), 
           d_size(0)
    {
	DBCA_Support::common_construct<T>(d_ptr, d_size, n, false);
    }

    //! Create an array of size \p n initialized to \p value.
    DBC_Array(const size_type n, const value_type& value)
        : d_ptr(0), 
          d_size(0)
    {
	DBCA_Support::common_construct<T>(d_ptr, d_size, n, true, value);
    }

    //! Create an array from the iterators [first, last)
    template<class InputIterator>
    DBC_Array(InputIterator first, InputIterator last)
        : d_ptr(0), 
          d_size(0)
    {
	// This dispatcher mechanism is used to differentiate between
	// InputIterator being an integral type or not.  If it is, it cannot
	// be an iterator, and so, should generate an array of length "first"
	// and initialized to value "last".
	using namespace DBCA_Support;
	initialize_dispatcher<T,InputIterator,
	    std::numeric_limits<InputIterator>::is_integer>::
	    do_it(d_ptr, d_size, first, last);
    }
    
    //! (Deep) Copy constructor
    DBC_Array(const DBC_Array<T>& rhs);

    //! Destructor
    ~DBC_Array()   
    { 
	if(d_ptr) delete[] d_ptr; 
    }

    //! (deep) copy
    DBC_Array<T>& operator=(const DBC_Array<T>& rhs);


    /*! 
      This has the same initialization behavior as DBC_Array(size_t), but,
      beyond that, there are \b no guarantees about the contents of the new
      array. 
    */
    void resize(const size_type n)
    {
	DBCA_Support::common_resize<T>(d_ptr, d_size, n, false);
    }

    //! Make the array of size \p n, initialized to \p value
    void assign(const size_type n, const_reference value)
    {
	DBCA_Support::common_resize<T>(d_ptr, d_size, n, true, value);
    }

    //! Change the array to be the same as the iterator range [first,last)
    template<class InputIterator>
    void assign(InputIterator first, InputIterator last)
    {
	// This dispatcher mechanism is used to differentiate between
	// InputIterator being an integral type or not.  If it is, it cannot
	// be an iterator, and so, should generate an array of length "first"
	// and initialized to value "last".
	using namespace DBCA_Support;
	assign_dispatcher<T,InputIterator,
	    std::numeric_limits<InputIterator>::is_integer>::
	    do_it(d_ptr, d_size, first, last);
    }

    //! Swap the contents of two arrays
    void swap(DBC_Array<T>& rhs);

    //! Release the memory (setting size to zero)
    void clear();

    //! Return the number of elements of type T that have been allocated
    size_type size() const  
    { 
	return d_size; 
    }

    //! Returns \c true if the array is empty
    bool empty() const 
    { 
	return d_size == 0; 
    }

    //! Returns a reference to the i'th element
    reference operator[](const size_t i) 
    {
	Require(i < d_size);
	return d_ptr[i];
    }

    //! Returns a const reference to the i'th element
    const_reference operator[](const size_t i) const 
    {
	Require(i < d_size);
	return d_ptr[i];
    }


    //! Returns a COPY of the i'th element
    value_type operator()(const size_t i) const
    {
	Require(i < d_size);
	return d_ptr[i];
    }

    //! Returns the begin iterator (T*)
    iterator begin()  
    { 
	return d_ptr; 
    }

    //! Returns the begin const iterator (T const *)
    const_iterator begin() const  
    { 
	return d_ptr; 
    }

    //! Returns the end iterator (T*): illegal to access contents here
    iterator end()  
    { 
	return d_ptr + d_size; 
    }

    //! Returns the end const iterator (T*): illegal to access contents here
    const_iterator end() const  
    { 
	return d_ptr + d_size; 
    }
    
    //! Returns a reference to the first element
    reference front() 
    {
	Require(d_size);
	return *(d_ptr);
    }

    //! Returns a const reference to the first element
    const_reference front() const
    {
	Require(d_size);
	return *(d_ptr);
    }

    //! Returns a reference to the last element
    reference back() 
    {
	Require(d_size);
	return *(d_ptr + d_size - 1);
    }

    //! Returns a const reference to the last element
    const_reference back() const
    {
	Require(d_size);
	return *(d_ptr + d_size - 1);
    }

    //! Element-wise equality
    bool operator==(const DBC_Array<T>& rhs) const;

    //! Element-wise less-than
    bool operator<(const DBC_Array<T>& rhs) const;

  private:
    //! Pointer to the real data
    T* d_ptr;
    //! Number of T allocated.
    size_type d_size;

};

//---------------------------------------------------------------------------//
//! Convenience output function
template<class T> std::ostream&
operator<<(std::ostream& os, const DBC_Array<T>& rhs);


//---------------------------------------------------------------------------//
// FREE COMPARISON FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//! Element-wise lhs != rhs
template<class T> inline bool
operator!=(const DBC_Array<T>& lhs, const DBC_Array<T>& rhs)
{
    return !(lhs == rhs);
}


//---------------------------------------------------------------------------//
//! Element-wise lhs > rhs
template<class T> inline bool
operator>(const DBC_Array<T>& lhs, const DBC_Array<T>& rhs)
{
    return rhs < lhs;
}


//---------------------------------------------------------------------------//
//! Element-wise lhs <= rhs
template<class T> inline bool
operator<=(const DBC_Array<T>& lhs, const DBC_Array<T>& rhs)
{
    return !(rhs < lhs);
}


//---------------------------------------------------------------------------//
//! Element-wise lhs >= rhs
template<class T> inline bool
operator>=(const DBC_Array<T>& lhs, const DBC_Array<T>& rhs)
{
    return !(lhs < rhs);
}


//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template<class T> 
DBC_Array<T>::DBC_Array(const DBC_Array<T>& rhs)
     : d_ptr(0), 
       d_size(0)
{
    DBCA_Support::common_construct<T>(d_ptr, d_size, rhs.d_size, false);
    if(d_size)
    {
	std::copy(rhs.begin(), rhs.end(), d_ptr);
    }
}


//---------------------------------------------------------------------------//
// OTHER PUBLIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
  If you are tempted to do something like memmove instead of copy, consider
  what happens when you copy things with non-trivial assignment operators.
  Many implementations of std::copy reduce to memmove if they can, anyway.
*/
template<class T> DBC_Array<T>&
DBC_Array<T>::operator=(const DBC_Array<T>& rhs) 
{
    if(&rhs != this) 
    {
	DBCA_Support::common_resize<T>(d_ptr, d_size, rhs.d_size, false);
	std::copy(rhs.begin(), rhs.end(), d_ptr);
    }
    return *this;
}


//---------------------------------------------------------------------------//
template<class T> void
DBC_Array<T>::swap(DBC_Array<T>& rhs) 
{
    if(&rhs != this)
    {
	std::swap(d_size, rhs.d_size);
	std::swap(d_ptr, rhs.d_ptr);
    }
}

//---------------------------------------------------------------------------//
template<class T> void 
DBC_Array<T>::clear() 
{
    if(d_size)
    {
	Check(d_ptr != 0);
	delete[] d_ptr;
	d_ptr = 0;
	d_size = 0;
    }
    Ensure(d_ptr == 0);
}

//---------------------------------------------------------------------------//
// COMPARISON FUNCTIONS
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
template<class T> bool
DBC_Array<T>::operator==(DBC_Array<T> const & rhs) const 
{
    if(rhs.d_size != d_size) return false;
    if(!d_size) return true;
    return std::equal(begin(), end(), rhs.d_ptr);
}


//---------------------------------------------------------------------------//
template<class T> bool
DBC_Array<T>::operator<(DBC_Array<T> const & rhs) const 
{
    return std::lexicographical_compare(begin(), end(),
					rhs.begin(), rhs.end());
}



//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
/*! 
  You could argue that this is superfluous, but, I just find it darn
  useful. 
*/
template<class T> std::ostream&
operator<<(std::ostream& os, 
	   const DBC_Array<T>& rhs)
{
    const size_t N = rhs.size();
    if(!N) return os;
    os << rhs[0];
    for(size_t i = 1; i < N; ++i)
	os << ' ' << rhs[i];
    return os;
}

}


#endif
