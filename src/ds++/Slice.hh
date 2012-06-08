//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   container/Slice.hh
 * \author Kent Budge
 * \date   Thu Jul  8 08:06:53 2004
 * \brief  Definition of Slice template
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef container_Slice_hh
#define container_Slice_hh

#include <iterator>
#include "Assert.hh"

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class Slice
 * \brief A reference to a Slice from a random access sequence.
 *
 * This class has reference semantics to the underlying random access
 * sequence, and imposes a subset on the sequence based on a stride.  It can
 * also be used to give standard container semantics to a sequence, which is
 * particularly useful with builtin types like double* that have iterator
 * semantics. 
 *
 * \arg \a Ran A random access iterator type, such as double* or
 * std::vector<double>::iterator.
 */
//===========================================================================//

template<class Ran>
class Slice 
{
  private:
    
    typedef typename std::iterator_traits<Ran> traits;

  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef typename traits::value_type      value_type;
    // no allocator type
    typedef unsigned size_type;
    typedef typename traits::difference_type difference_type;

    class iterator
    {
      public:
	typedef typename traits::iterator_category iterator_category;
	typedef typename traits::value_type        value_type;
	typedef typename traits::difference_type   difference_type;
	typedef typename traits::pointer           pointer;
	typedef typename traits::reference         reference;

	iterator &operator++()
	{
	    offset += stride;
	    return *this;
	}

	iterator operator++(int)
	{
	    offset += stride;
	    return iterator(first, offset-stride, stride);
	}

	reference operator*() const
	{
	    return first[offset];
	}

	reference operator[](difference_type i) const
	{
	    return first[offset+i];
	}

	iterator operator+(difference_type i) const
	{
	    return iterator(first, offset+i*stride, stride);
	}

	difference_type operator-(iterator i) const
	{
	    Require((first-i.first)%stride==0);
	    Require(stride==i.stride);

	    return ((first-i.first)+offset-i.offset)/stride;
	}

      private:
	friend class Slice;
        friend class const_iterator;

	iterator(Ran const first,
		 difference_type const offset,
		 unsigned const stride)
	    : first(first),
	      offset(offset),
	      stride(stride)
	{
	    Require(stride>0);
	}

	friend bool operator<(iterator a, iterator b)
	{
	    return a.first-b.first<b.offset-a.offset;
	}

	friend bool operator!=(iterator a, iterator b)
	{
	    return a.first-b.first!=b.offset-a.offset;
	}

	Ran first;
	difference_type offset;
	unsigned stride;
    };

    class const_iterator
    {
      public:
	typedef typename traits::iterator_category iterator_category;
	typedef typename traits::value_type        value_type;
	typedef typename traits::difference_type   difference_type;
	typedef typename traits::pointer           pointer;
	typedef typename traits::reference         reference;

	const_iterator &operator++()
	{
	    offset += stride;
	    return *this;
	}

	const_iterator operator++(int)
	{
	    offset += stride;
	    return const_iterator(first, offset-stride, stride);
	}

	value_type const & operator*() const
	{
	    return first[offset];
	}

	value_type const & operator[](difference_type i) const
	{
	    return first[offset+i];
	}

	const_iterator operator+(difference_type i) const
	{
 	    return const_iterator(first, offset+i*stride, stride);
	}

	difference_type operator-(const_iterator i) const
	{
	    Require((first-i.first)%stride==0);
	    Require(stride==i.stride);

	    return ((first-i.first)+offset-i.offset)/stride;
	}

	const_iterator(iterator const &i)
	    : first(i.first),
	      offset(i.offset),
	      stride(i.stride)
	{
	    Require(stride>0);
	}

      private:
	friend class Slice;

	const_iterator(Ran const first,
		       difference_type const offset,
		       unsigned const stride)
	    : first(first),
	      offset(offset),
	      stride(stride)
	{
	    Require(stride>0);
	}

	friend bool operator<(const_iterator a, const_iterator b)
	{
	    return a.first-b.first<b.offset-a.offset;
	}

	friend bool operator!=(const_iterator a, const_iterator b)
	{
	    return a.first-b.first!=b.offset-a.offset;
	}

	Ran first;
	difference_type offset;
	unsigned stride;
    };

    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    typedef typename traits::pointer   pointer;
    typedef typename traits::reference reference;
    typedef value_type const *         const_pointer;
    typedef value_type const &         const_reference;

    // CREATORS

    //-----------------------------------------------------------------------//
    /*! 
     * \brief Construct from a sequence. 
     * 
     * The constructed Slice has reference semantics to the sequence.  That
     * is, modifications to elements of the Slice are actually modifications
     * to the elements of the underlying sequence, and the Slice becomes
     * invalid if the underlying sequence becomes invalid.  For example, a
     * sequence based on a vector becomes invalid if the vector is resized,
     * and so will a Slice based on that sequence.
     *
     * A stride is applied to all indexing, so that an index \c i applied to
     * the Slice references the element whose index is \c stride*i in the
     * underlying sequence.  No offset is necessary, since the same result is
     * easily obtained by modifying the starting iterator of the underlying
     * sequence. Iterators behave in a comparable manner, so
     * that \c begin() of the Slice points to \c begin()+offset of the
     * underlying container, and the increment operator actually increments
     * by \c stride.
     *
     * \param first Iterator to the beginning of a sequence.
     * \param length Length of the constructed Slice.
     * \param stride Stride to apply to the sequence.
     */
    Slice(Ran const first,
	  unsigned const length,
	  unsigned const stride = 1)
	: first(first), length(length), stride(stride)
    {
	Require(stride>0);
    }

    // MANIPULATORS

    // ACCESSORS

    iterator begin()
    {
	return iterator(first, 0, stride);
    }
    const_iterator begin() const
    {
	return const_iterator(first, 0, stride);
    }
    iterator end()
    {
	return iterator(first, length*stride, stride);
    }
    const_iterator end() const
    {
	return const_iterator(first, length*stride, stride);
    }

    reverse_iterator rbegin();
    const_reverse_iterator rbegin() const;
    reverse_iterator rend();
    const_reverse_iterator rend() const;

    reference operator[](size_type n);
    const_reference operator[](size_type n) const;

    reference at(size_type n);
    const_reference at(size_type n) const;

    reference front();
    const_reference front() const;
    reference back();
    const_reference back() const;

    size_type size() const { return length; }
    bool empty() const { return size()==0; }
    size_type max_size() const;

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA

    Ran first;
    unsigned length;
    unsigned stride;
};

//---------------------------------------------------------------------------//
/*! 
 * \brief Return a Slice.
 *
 * This is a convenient interface that avoids the necessity for writing a lot
 * of type expressions in user code.
 * 
 * \arg \a Rab A random access iterator.
 *
 * \param first Start of a sequence
 * \param last Length of the Slice.
 * \param stride Stride into the container
 *
 * \return The desired Slice.
 */

template<class Ran>
inline
Slice<Ran> slice(Ran const first,
		 unsigned const length,
		 unsigned const stride = 1)
{
    Require(stride>0);

    return Slice<Ran>(first, length, stride);
}

} // end namespace rtt_dsxx

#include "Slice.t.hh"

#endif // container_Slice_hh

//---------------------------------------------------------------------------//
//              end of container/Slice.hh
//---------------------------------------------------------------------------//
