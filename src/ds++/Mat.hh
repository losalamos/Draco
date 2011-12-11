//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Mat.hh
 * \author Geoffrey Furnish
 * \date   Fri Jan 24 15:48:31 1997
 * \brief  Mat class definitions.
 * \note   Copyright (C) 2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_ds_Mat_hh
#define rtt_ds_Mat_hh

#include "Assert.hh"
#include "Bounds.hh"
#include "destroy.hh"
#include <algorithm>
#include <memory>

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class Mat1
 * \brief A 1-d container.
 *
 * This class is intended to be a replacement for the STL vector<T> class.
 * The reasons for wanting a replacement for vector<T> are primarily
 * threefold: 
 * -# Integration with the DS++ assertion model. 
 * -# Ability to provide more sophisticated allocators.
 * -# provide a 1-d analog to Mat2 and Mat3.
 * .
 */
/*!
 * \example ds++/test/tstMat1RA.cc
 * 
 * Test of Mat1.
 */
//===========================================================================//

template< typename T,
	  typename Allocator = typename std::allocator<T> >
//	  class Allocator = typename alloc_traits<T>::Default_Allocator >
class Mat1
{
  private:
    int xmin, xlen;
    bool may_free_space;

    int xmax() const { return xmin + xlen - 1; }

    // index() is used for indexing, so is checked.  offset() is used for
    // adjusting the memory pointer, which is logically distinct from
    // indexing, so is not checked.  Note in particular, that a zero size
    // matrix will have no valid index, but we still need to be able to
    // compute the offsets.

    int index( int i ) const
    {
	Assert( i >= xmin );
	Assert( i < xmin + xlen );
	return offset(i);
    }
    int offset( int i ) const { return i; }

    void detach()
    {
	if (may_free_space)
        {
	    rtt_dsxx::Destroy( begin(), end() );
	    alloc.deallocate( v + offset(xmin), size() );
	}
    }

  protected:
    Allocator alloc;
    T *v;

  public:
    typedef       T  value_type;
    typedef       T& reference;
    typedef const T& const_reference;
    typedef       T* pointer;
    typedef const T* const_pointer;
    typedef typename Allocator::difference_type difference_type;
    typedef typename Allocator::size_type       size_type;
    typedef typename Allocator::pointer         iterator;
    typedef typename Allocator::const_pointer   const_iterator;
    typedef typename std::reverse_iterator<iterator>       reverse_iterator;
    typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;

    // Accessors

    T&       operator[]( int i )       { return v[ index(i) ]; }
    const T& operator[]( int i ) const { return v[ index(i) ]; }

    T&       operator()( int i )       { return v[ index(i) ]; }
    const T& operator()( int i ) const { return v[ index(i) ]; }

    iterator       begin()       { return v + offset(xmin); }
    const_iterator begin() const { return v + offset(xmin); }

    iterator       end()         { return v + offset(xmax()) + 1; }
    const_iterator end() const   { return v + offset(xmax()) + 1; }

    reverse_iterator       rbegin()       { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(end()); }

    reverse_iterator       rend()         { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(begin()); }

    int get_xmin() const { return xmin; }
    int get_xlen() const { return xlen; }

    size_type size() const { return nx(); }
    size_type max_size () const { return alloc.max_size(); }
    bool empty() const { return (this->size() == 0); }

    // For backward compatibility.
    int nx() const { return xlen; }

    // Constructors
    Mat1(void)
	: xmin(0),
          xlen(0),
	  may_free_space(false),
          alloc(),
          v(0)
    {/*empty*/}

    explicit Mat1( int xmax_, const T& t = T() )
	: xmin(0),
          xlen(xmax_),
	  may_free_space(true),
	  alloc(),
          v( alloc.allocate( size() ) - offset(xmin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat1( T *vv, int xmax_ )
	: xmin(0),
          xlen(xmax_),
	  may_free_space(false),
          alloc(),
          v(vv)
    {/*empty*/}

    explicit Mat1( const Bounds& bx, const T& t = T() )
	: xmin( bx.min() ),
          xlen( bx.len() ),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

//     Mat1( T *vv, const Bounds& bx )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  may_free_space(false),
// 	  v( vv - offset(xmin) )
//     {}

    Mat1( const Mat1<T>& m )
	: xmin(m.xmin),
          xlen(m.xlen),
	  may_free_space(true),
	  alloc(),
          v( alloc.allocate( size() ) - offset(xmin) )
    {
	std::uninitialized_copy( m.begin(), m.end(), begin() );
    }

    // Destructor

    ~Mat1()
    {
	detach();
    }

    // Assignment operators

    Mat1& operator=( const T& t )
    {
	std::fill( begin(), end(), t );
	return *this;
    }

    Mat1& operator=( const Mat1& m )
    {
	if (this == &m) return *this;

	if ( m.xmin != xmin || m.xlen != xlen ) {
	    detach();
	    xmin = m.xmin;
	    xlen = m.xlen;
	    v = alloc.allocate( size() ) - offset(xmin);
	    std::uninitialized_copy( m.begin(), m.end(), begin() );
	    may_free_space = true;
	}
	else {
	    if (v)
		std::copy( m.begin(), m.end(), begin() );
	}

	return *this;
    }

    void swap ( Mat1& m )
    {
        int itemp;
        bool btemp;
        Allocator atemp;
        T* ptemp;

        itemp = xmin;
        xmin = m.xmin;
        m.xmin = itemp;

        itemp = xlen;
        xlen = m.xlen;
        m.xlen = itemp;

        btemp = may_free_space;
        may_free_space = m.may_free_space;
        m.may_free_space = btemp;

        atemp = alloc;
        alloc = m.alloc;
        m.alloc = atemp;

        ptemp = v;
        v = m.v;
        m.v = ptemp;
    }

    // Boolean operators

    bool operator==( const Mat1& m ) const
    {
        if (this == &m)
            return true;

        if ( m.size() != this->size() )
            return false;

        const_iterator miter = m.begin();
        for (const_iterator iter = this->begin(); iter != this->end(); ++iter)
        {
            if (*iter != *miter) return false;
            ++miter;
        }

        return true;
    }

    bool operator!=( const Mat1& m ) const
    {
        return !(*this == m);
    }

    bool operator<( const Mat1& m ) const
    {
        return std::lexicographical_compare(this->begin(), this->end(), 
					    m.begin(), m.end());
    }

    bool operator>( const Mat1& m ) const
    {
        return (m < *this);
    }

    bool operator<=( const Mat1& m ) const
    {
        return !(m < *this);
    }

    bool operator>=( const Mat1& m ) const
    {
        return !(*this < m);
    }

    // Mathematical support

    template<class X> Mat1& operator+=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ += x;
	return *this;
    }
    Mat1& operator+=( const Mat1<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ += *j++;
	return *this;
    }

    template<class X> Mat1& operator-=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ -= x;
	return *this;
    }
    Mat1& operator-=( const Mat1<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ -= *j++;
	return *this;
    }

    template<class X> Mat1& operator*=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ *= x;
	return *this;
    }
    Mat1& operator*=( const Mat1<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ *= *j++;
	return *this;
    }

    template<class X> Mat1& operator/=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ /= x;
	return *this;
    }
    Mat1& operator/=( const Mat1<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ /= *j++;
	return *this;
    }

    // Utility support

    void assert_conformality( const Mat1<T>& Remember(m) ) const
    {
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
    }

    void redim( int nxmax, const T& t = T() )
    {
	// This one only works right if xmin == 0.
	Assert( xmin == 0 );
	if (v && !may_free_space) {
	    // User thinks he wants to expand the aliased region.
	    xlen = nxmax;
	    return;
	}
	detach();
	xlen = nxmax;
	v = alloc.allocate( size() ) - offset(xmin);
	std::uninitialized_fill( begin(), end(), t );
	may_free_space = true;
    }

    void redim( const Bounds& bx, const T& t = T() )
    {
	if (v && !may_free_space) {
	    // Respecify the aliased region.
	    v += offset(xmin);
	    xmin = bx.min(); xlen = bx.len();
	    v -= offset(xmin);
	    return;
	}
	detach();
	xmin = bx.min(); xlen = bx.len();
	v = alloc.allocate( size() ) - offset(xmin);
	std::uninitialized_fill( begin(), end(), t );
	may_free_space = true;
    }

    // Check to see if this Mat1<T> is of size x.

    bool conformal( int x ) const
    {
	return xmin == 0 && x == xlen;
    }

    // Obtain dimension of this Mat1<T>.

    void elements( int& nx_ ) const { nx_ = nx(); }
};

//===========================================================================//
/*!
 * \class Mat2 
 * \brief A 2-d container.
 *
 * Mat2 provides a container which can be indexed using two indices.  The STL
 * provides nothing quite like this, but it is essential for many
 * mathematical purposes to have a 2-d container.
 */
/*!
 * \example ds++/test/tstMat2RA.cc
 * 
 * Test of Mat2.
 */
//===========================================================================//

template< typename T,
	  typename Allocator = typename std::allocator<T> >          
          //class Allocator = Simple_Allocator<T> >
class Mat2
{
  private:
    int xmin, xlen, ymin, ylen;
    bool may_free_space;

    int xmax() const { return xmin + xlen - 1; }
    int ymax() const { return ymin + ylen - 1; }

    // index() is used for indexing, so is checked.  offset() is used for
    // adjusting the memory pointer, which is logically distinct from
    // indexing, so is not checked.  Note in particular, that a zero size
    // matrix will have no valid index, but we still need to be able to
    // compute the offsets.

    // Compute the offset into the data array, of the i,j th element.
    int index( int i, int j ) const
    {
	Assert( i >= xmin );
	Assert( i < xmin + xlen );
	Assert( j >= ymin );
	Assert( j < ymin + ylen );

	return offset(i,j);
    }
    int offset( int i, int j ) const { return xlen * j + i; }

    // Make sure a bare integer index is within the appropriate range.
    void check( int Remember(i) ) const
    {
	Assert( i >= offset( xmin, ymin ) );
	Assert( i <= offset( xmax(), ymax() ) );
    }

    void detach()
    {
	if (may_free_space) {
	    rtt_dsxx::Destroy( begin(), end() );
	    alloc.deallocate( v + offset(xmin,ymin), size() );
	}
    }

  protected:
    Allocator alloc;
    T *v;

  public:
    typedef       T  value_type;
    typedef       T& reference;
    typedef const T& const_reference;
    typedef       T* pointer;
    typedef const T* const_pointer;
    typedef typename Allocator::difference_type difference_type;
    typedef typename Allocator::size_type       size_type;
    typedef typename Allocator::pointer         iterator;
    typedef typename Allocator::const_pointer   const_iterator;
    typedef typename std::reverse_iterator<iterator>       reverse_iterator;
    typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;

    // Accessors

    T&       operator()( int i, int j )       { return v[ index(i,j) ]; }
    const T& operator()( int i, int j ) const { return v[ index(i,j) ]; }

    T& operator[]( int i ) { check(i); return v[i]; }
    const T& operator[]( int i ) const { check(i); return v[i]; }

    iterator       begin()       { return v + offset(xmin,ymin); }
    const_iterator begin() const { return v + offset(xmin,ymin); }

    iterator       end()         { return v + offset(xmax(),ymax()) + 1; }
    const_iterator end() const   { return v + offset(xmax(),ymax()) + 1; }

    reverse_iterator       rbegin()       { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(end()); }

    reverse_iterator       rend()         { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(begin()); }

    int get_xmin() const { return xmin; }
    int get_xlen() const { return xlen; }
    int get_ymin() const { return ymin; }
    int get_ylen() const { return ylen; }

    size_type size() const { return nx() * ny(); }
    size_type max_size () const { return alloc.max_size(); }
    bool empty() const { return (this->size() == 0); }

    // For backward compatibility.
    int nx() const { return xlen; }
    int ny() const { return ylen; }

    // Constructors

    Mat2(void)
	: xmin(0),
          xlen(0),
          ymin(0),
          ylen(0),
	  may_free_space(false),
          alloc(),
          v(0)
    {/*empty*/}

    Mat2( int xmax_, int ymax_, const T& t = T() )
	: xmin(0),
          xlen(xmax_),
	  ymin(0),
          ylen(ymax_),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat2( T *vv, int xmax_, int ymax_ )
	: xmin(0),
          xlen(xmax_),
	  ymin(0),
          ylen(ymax_),
	  may_free_space(false),
          alloc(),
          v(vv)
    {/*empty*/}

    Mat2( const Bounds& bx, const Bounds& by, const T& t = T() )
	: xmin( bx.min() ),
          xlen( bx.len() ),
	  ymin( by.min() ),
          ylen( by.len() ),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

// [kt - 3/6/2006] No unit tests for this function.  Since it is not used, we
// don't need to define it.
//     Mat2( T *vv, const Bounds& bx, const Bounds& by )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  ymin( by.min() ), ylen( by.len() ),
// 	  may_free_space(false),
// 	  v( vv - offset(xmin,ymin) )
//     {}

    Mat2( const Mat2<T>& m )
	: xmin(m.xmin),
          xlen(m.xlen),
	  ymin(m.ymin),
          ylen(m.ylen),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin) )
    {
	std::uninitialized_copy( m.begin(), m.end(), begin() );
    }

    // Destructor

    ~Mat2()
    {
	detach();
    }

    // Assignment operators

    Mat2& operator=( const T& t )
    {
	std::fill( begin(), end(), t );
	return *this;
    }

    Mat2& operator=( const Mat2& m )
    {
	if (this == &m) return *this;

	if ( m.xmin != xmin || m.xlen != xlen ||
	     m.ymin != ymin || m.ylen != ylen ) {
	    detach();
	    xmin = m.xmin;
	    xlen = m.xlen;
	    ymin = m.ymin;
	    ylen = m.ylen;
	    v = alloc.allocate( size() ) - offset(xmin,ymin);
	    std::uninitialized_copy( m.begin(), m.end(), begin() );
	    may_free_space = true;
	}
	else {
	    if (v)
		std::copy( m.begin(), m.end(), begin() );
	}

	return *this;
    }

    void swap ( Mat2& m )
    {
        int itemp;
        bool btemp;
        Allocator atemp;
        T* ptemp;

        itemp = xmin;
        xmin = m.xmin;
        m.xmin = itemp;

        itemp = xlen;
        xlen = m.xlen;
        m.xlen = itemp;

        itemp = ymin;
        ymin = m.ymin;
        m.ymin = itemp;

        itemp = ylen;
        ylen = m.ylen;
        m.ylen = itemp;

        btemp = may_free_space;
        may_free_space = m.may_free_space;
        m.may_free_space = btemp;

        atemp = alloc;
        alloc = m.alloc;
        m.alloc = atemp;

        ptemp = v;
        v = m.v;
        m.v = ptemp;
    }

    // Boolean operators

    bool operator==( const Mat2& m ) const
    {
        if (this == &m)
            return true;

        if ( m.size() != this->size() )
            return false;

        const_iterator miter = m.begin();
        for (const_iterator iter = this->begin(); iter != this->end(); ++iter)
        {
            if (*iter != *miter) return false;
            ++miter;
        }

        return true;
    }

    bool operator!=( const Mat2& m ) const
    {
        return !(*this == m);
    }

    bool operator<( const Mat2& m ) const
    {
        return std::lexicographical_compare(this->begin(), this->end(), 
					    m.begin(), m.end());
    }

    bool operator>( const Mat2& m ) const
    {
        return (m < *this);
    }

    bool operator<=( const Mat2& m ) const
    {
        return !(m < *this);
    }

    bool operator>=( const Mat2& m ) const
    {
        return !(*this < m);
    }

    // Mathematical support

    template<class X> Mat2& operator+=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ += x;
	return *this;
    }
    Mat2& operator+=( const Mat2<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ += *j++;
	return *this;
    }

    template<class X> Mat2& operator-=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ -= x;
	return *this;
    }
    Mat2& operator-=( const Mat2<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ -= *j++;
	return *this;
    }

    template<class X> Mat2& operator*=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ *= x;
	return *this;
    }
    Mat2& operator*=( const Mat2<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ *= *j++;
	return *this;
    }

    template<class X> Mat2& operator/=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ /= x;
	return *this;
    }
    Mat2& operator/=( const Mat2<T>& m )
    {
	assert_conformality( m );
	iterator i = begin();
	const_iterator j = m.begin();
	while( i != end() ) *i++ /= *j++;
	return *this;
    }

    // Utility support

    void assert_conformality( const Mat2<T>& Remember(m) ) const
    {
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
	Assert( ymin == m.ymin );
	Assert( ylen == m.ylen );
    }

    void redim( int nxmax, int nymax, const T& t = T() )
    {
	// This one only works right if xmin == 0 and ymin == 0.
	Assert( xmin == 0 );
	Assert( ymin == 0 );
	if (v && !may_free_space) {
	    // User thinks he wants to expand the aliased region.
	    xlen = nxmax;
	    ylen = nymax;
	    return;
	}
	detach();
	xlen = nxmax;
	ylen = nymax;
	v = alloc.allocate( size() ) - offset(0,0);
	std::uninitialized_fill( begin(), end(), t );
	may_free_space = true;
    }

    void redim( const Bounds& bx, const Bounds& by, const T& t = T() )
    {
	if (v && !may_free_space) {
	    // Respecify the aliased region.
	    v += offset(xmin,ymin);
	    xmin = bx.min(); xlen = bx.len();
	    ymin = by.min(); ylen = by.len();
	    v -= offset(xmin,ymin);
	    return;
	}
	detach();
	xmin = bx.min(); xlen = bx.len();
	ymin = by.min(); ylen = by.len();
	v = alloc.allocate( size() ) - offset(xmin,ymin);
	std::uninitialized_fill( begin(), end(), t );
	may_free_space = true;
    }

    // Check to see if this Mat2<T> is of size x by y.

    bool conformal( int x, int y ) const
    {
	return xmin == 0 && x == xlen &&
	    ymin == 0 && y == ylen;
    }

    // Obtain dimensions of this Mat2<T>.

    void elements( int& nx_, int& ny_ ) const { nx_ = nx(); ny_ = ny(); }
};

//===========================================================================//
/*!
 * \class Mat3
 * \brief A 3-d container.
 *
 * Mat3 is a container which supports three indices.  It is otherwise similar
 * to Mat1 and Mat2.
 */
/*!
 * \example ds++/test/tstMat3RA.cc
 * 
 * Test of Mat3.
 */
//===========================================================================//

template< typename T,
	  typename Allocator = typename std::allocator<T> >
class Mat3
{
  private:
    int xmin, xlen, ymin, ylen, zmin, zlen;
    bool may_free_space;

    int xmax() const { return xmin + xlen - 1; }
    int ymax() const { return ymin + ylen - 1; }
    int zmax() const { return zmin + zlen - 1; }

    // index() is used for indexing, so is checked.  offset() is used for
    // adjusting the memory pointer, which is logically distinct from
    // indexing, so is not checked.  Note in particular, that a zero size
    // matrix will have no valid index, but we still need to be able to
    // compute the offsets.

    // Compute the offset into the data array, of the i,j th element.
//     int index( int i, int j, int k ) const
//     {
// 	Assert( i >= xmin );
// 	Assert( i < xmin + xlen );
// 	Assert( j >= ymin );
// 	Assert( j < ymin + ylen );
// 	Assert( k >= zmin );
// 	Assert( k < zmin + zlen );

// 	return offset(i,j,k);
//     }
    int offset( int i, int j, int k ) const { return xlen*(k*ylen+j)+i; }

    // Make sure a bare integer index is within the appropriate range.
    void check( int Remember(i) ) const
    {
	Assert( i >= offset( xmin, ymin, zmin ) );
	Assert( i <= offset( xmax(), ymax(), zmax() ) );
    }

    void detach()
    {
	if (may_free_space) {
	    rtt_dsxx::Destroy( begin(), end() );
	    alloc.deallocate( v + offset(xmin,ymin,zmin), size() );
	}
    }

  protected:
    Allocator alloc;
    T *v;

  public:
    typedef       T  value_type;
    typedef       T& reference;
    typedef const T& const_reference;
    typedef       T* pointer;
    typedef const T* const_pointer;
    typedef typename Allocator::difference_type difference_type;
    typedef typename Allocator::size_type       size_type;
    typedef typename Allocator::pointer         iterator;
    typedef typename Allocator::const_pointer   const_iterator;
    typedef typename std::reverse_iterator<iterator>       reverse_iterator;
    typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;
    
    // Accessors

//     T&       operator()( int i, int j, int k )       { return v[ index(i,j,k) ]; }
//     const T& operator()( int i, int j, int k ) const { return v[ index(i,j,k) ]; }

    T& operator[]( int i ) { check(i); return v[i]; }
    const T& operator[]( int i ) const { check(i); return v[i]; }

    iterator       begin()       { return v + offset(xmin,ymin,zmin); }
    const_iterator begin() const { return v + offset(xmin,ymin,zmin); }

    iterator       end()         { return v + offset(xmax(),ymax(),zmax()) + 1; }
    const_iterator end() const   { return v + offset(xmax(),ymax(),zmax()) + 1; }

    reverse_iterator       rbegin()       { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(end()); }

    reverse_iterator       rend()         { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(begin()); }

//     int get_xmin() const { return xmin; }
//     int get_xlen() const { return xlen; }
//     int get_ymin() const { return ymin; }
//     int get_ylen() const { return ylen; }
//     int get_zmin() const { return zmin; }
//     int get_zlen() const { return zlen; }

    size_type size() const { return nx() * ny() * nz(); }
    size_type max_size () const { return alloc.max_size(); }
    bool empty() const { return (this->size() == 0); }

    // For backward compatibility.
    int nx() const { return xlen; }
    int ny() const { return ylen; }
    int nz() const { return zlen; }

    // Constructors

//     Mat3()
// 	: xmin(0), xlen(0), ymin(0), ylen(0), zmin(0), zlen(0),
// 	  may_free_space(false), v(0)
//     {}

    Mat3( int xmax_, int ymax_, int zmax_, const T& t = T() )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  zmin(0), zlen(zmax_),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat3( T *vv, int xmax_, int ymax_, int zmax_ )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  zmin(0), zlen(zmax_),
	  may_free_space(false),
          alloc(),
          v(vv)
    {/*empty*/}

//     Mat3( const Bounds& bx, const Bounds& by,
// 	  const Bounds& bz, const T& t = T() )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  ymin( by.min() ), ylen( by.len() ),
// 	  zmin( bz.min() ), zlen( bz.len() ),
// 	  may_free_space(true),
// 	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin) )
//     {
// 	std::uninitialized_fill( begin(), end(), t );
//     }

//     Mat3( T *vv, const Bounds& bx, const Bounds& by, const Bounds& bz )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  ymin( by.min() ), ylen( by.len() ),
// 	  zmin( bz.min() ), zlen( bz.len() ),
// 	  may_free_space(false),
// 	  v( vv - offset(xmin,ymin,zmin) )
//     {}

    Mat3( const Mat3<T>& m )
	: xmin(m.xmin), xlen(m.xlen),
	  ymin(m.ymin), ylen(m.ylen),
	  zmin(m.zmin), zlen(m.zlen),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin) )
    {
	std::uninitialized_copy( m.begin(), m.end(), begin() );
    }

    // Destructor

    ~Mat3()
    {
	detach();
    }

    // Assignment operators

//     Mat3& operator=( const T& t )
//     {
// 	std::fill( begin(), end(), t );
// 	return *this;
//     }

    Mat3& operator=( const Mat3& m )
    {
	if (this == &m) return *this;

	if ( m.xmin != xmin || m.xlen != xlen ||
	     m.ymin != ymin || m.ylen != ylen ||
	     m.zmin != zmin || m.zlen != zlen ) {
	    detach();
	    xmin = m.xmin;
	    xlen = m.xlen;
	    ymin = m.ymin;
	    ylen = m.ylen;
	    zmin = m.zmin;
	    zlen = m.zlen;
	    v = alloc.allocate( size() ) - offset(xmin,ymin,zmin);
	    std::uninitialized_copy( m.begin(), m.end(), begin() );
	    may_free_space = true;
	}
	else {
	    if (v)
		std::copy( m.begin(), m.end(), begin() );
	}

	return *this;
    }

    void swap ( Mat3& m )
    {
        int itemp;
        bool btemp;
        Allocator atemp;
        T* ptemp;

        itemp = xmin;
        xmin = m.xmin;
        m.xmin = itemp;

        itemp = xlen;
        xlen = m.xlen;
        m.xlen = itemp;

        itemp = ymin;
        ymin = m.ymin;
        m.ymin = itemp;

        itemp = ylen;
        ylen = m.ylen;
        m.ylen = itemp;

        itemp = zmin;
        zmin = m.zmin;
        m.zmin = itemp;

        itemp = zlen;
        zlen = m.zlen;
        m.zlen = itemp;

        btemp = may_free_space;
        may_free_space = m.may_free_space;
        m.may_free_space = btemp;

        atemp = alloc;
        alloc = m.alloc;
        m.alloc = atemp;

        ptemp = v;
        v = m.v;
        m.v = ptemp;
    }

    // Boolean operators

    bool operator==( const Mat3& m ) const
    {
        if (this == &m)
            return true;

        if ( m.size() != this->size() )
            return false;

        const_iterator miter = m.begin();
        for (const_iterator iter = this->begin(); iter != this->end(); ++iter)
        {
            if (*iter != *miter) return false;
            ++miter;
        }

        return true;
    }

    bool operator!=( const Mat3& m ) const
    {
        return !(*this == m);
    }

    bool operator<( const Mat3& m ) const
    {
        return std::lexicographical_compare(this->begin(), this->end(), 
					    m.begin(), m.end());
    }

    bool operator>( const Mat3& m ) const
    {
        return (m < *this);
    }

    bool operator<=( const Mat3& m ) const
    {
        return !(m < *this);
    }

    bool operator>=( const Mat3& m ) const
    {
        return !(*this < m);
    }

    // Mathematical support

//     template<class X> Mat3& operator+=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ += x;
// 	return *this;
//     }
//     Mat3& operator+=( const Mat3<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ += *j++;
// 	return *this;
//     }

//     template<class X> Mat3& operator-=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ -= x;
// 	return *this;
//     }
//     Mat3& operator-=( const Mat3<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ -= *j++;
// 	return *this;
//     }

//     template<class X> Mat3& operator*=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ *= x;
// 	return *this;
//     }
//     Mat3& operator*=( const Mat3<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ *= *j++;
// 	return *this;
//     }

//     template<class X> Mat3& operator/=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ /= x;
// 	return *this;
//     }
//     Mat3& operator/=( const Mat3<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ /= *j++;
// 	return *this;
//     }

    // Utility support

//     void assert_conformality( const Mat3<T>& m ) const
//     {
// 	Assert( xmin == m.xmin );
// 	Assert( xlen == m.xlen );
// 	Assert( ymin == m.ymin );
// 	Assert( ylen == m.ylen );
// 	Assert( zmin == m.zmin );
// 	Assert( zlen == m.zlen );
//     }

//     void redim( int nxmax, int nymax, int nzmax, const T& t = T() )
//     {
// 	// This one only works right if xmin == 0 and ymin == 0 and zmin == 0.
// 	Assert( xmin == 0 );
// 	Assert( ymin == 0 );
// 	Assert( zmin == 0 );
// 	if (v && !may_free_space) {
// 	    // User thinks he wants to expand the aliased region.
// 	    xlen = nxmax;
// 	    ylen = nymax;
// 	    zlen = nzmax;
// 	    return;
// 	}
// 	detach();
// 	xlen = nxmax;
// 	ylen = nymax;
// 	zlen = nzmax;
// 	v = alloc.allocate( size() ) - offset(0,0,0);
// 	std::uninitialized_fill( begin(), end(), t );
// 	may_free_space = true;
//     }

//     void redim( const Bounds& bx, const Bounds& by,
// 		const Bounds& bz, const T& t = T() )
//     {
// 	if (v && !may_free_space) {
// 	    // Respecify the aliased region.
// 	    v += offset(xmin,ymin,zmin);
// 	    xmin = bx.min(); xlen = bx.len();
// 	    ymin = by.min(); ylen = by.len();
// 	    zmin = bz.min(); zlen = bz.len();
// 	    v -= offset(xmin,ymin,zmin);
// 	    return;
// 	}
// 	detach();
// 	xmin = bx.min(); xlen = bx.len();
// 	ymin = by.min(); ylen = by.len();
// 	zmin = bz.min(); zlen = bz.len();
// 	v = alloc.allocate( size() ) - offset(xmin,ymin,zmin);
// 	std::uninitialized_fill( begin(), end(), t );
// 	may_free_space = true;
//     }

    // WARNING: This doesn't make a lot of sense anymore.

    // Check to see if this Mat3<T> is of size x by y by z.

//     bool conformal( int x, int y, int z ) const
//     {
// 	return x == xlen && y == ylen && z == zlen;
//     }

    // Obtain dimensions of this Mat3<T>.

//     void elements( int& nx, int& ny, int& nz ) const
//     {
// 	nx = xlen; ny = ylen; nz = zlen;
//     }
};

//===========================================================================//
/*!
 * \class Mat4
 * \brief A 4-d container.
 *
 * Mat3 is a container which supports four indices.  It is otherwise similar
 * to Mat1, Mat2 and Mat3.
 */
/*!
 * \example ds++/test/tstMat4RA.cc
 * 
 * Test of Mat4.
 */
//===========================================================================//

template< typename T,
	  typename Allocator = typename std::allocator<T> >
class Mat4
{
  private:
    int xmin, xlen, ymin, ylen, zmin, zlen, wmin, wlen;
    bool may_free_space;

    int xmax() const { return xmin + xlen - 1; }
    int ymax() const { return ymin + ylen - 1; }
    int zmax() const { return zmin + zlen - 1; }
    int wmax() const { return wmin + wlen - 1; }

    // index() is used for indexing, so is checked.  offset() is used for
    // adjusting the memory pointer, which is logically distinct from
    // indexing, so is not checked.  Note in particular, that a zero size
    // matrix will have no valid index, but we still need to be able to
    // compute the offsets.

    // Compute the offset into the data array, of the i,j,k,l th element.
//     int index( int i, int j, int k, int l ) const
//     {
// 	Assert( i >= xmin );
// 	Assert( i <  xmin + xlen );
// 	Assert( j >= ymin );
// 	Assert( j <  ymin + ylen );
// 	Assert( k >= zmin );
// 	Assert( k <  zmin + zlen );
// 	Assert( l >= wmin );
// 	Assert( l <  wmin + wlen );

// 	return offset(i,j,k,l);
//     }
    int offset( int i, int j, int k, int l ) const
    {
	return xlen*(ylen*(zlen*l+k)+j)+i;
    }

    // Make sure a bare integer index is within the appropriate range.
    void check( int Remember(i) ) const
    {
	Assert( i >= offset( xmin,   ymin,   zmin,   wmin   ) );
	Assert( i <= offset( xmax(), ymax(), zmax(), wmax() ) );
    }

    void detach()
    {
	if (may_free_space) {
	    rtt_dsxx::Destroy( begin(), end() );
	    alloc.deallocate( v + offset(xmin,ymin,zmin,wmin), size() );
	}
    }

  protected:
    Allocator alloc;
    T *v;

  public:
    typedef       T  value_type;
    typedef       T& reference;
    typedef const T& const_reference;
    typedef       T* pointer;
    typedef const T* const_pointer;
    typedef typename Allocator::difference_type difference_type;
    typedef typename Allocator::size_type       size_type;
    typedef typename Allocator::pointer         iterator;
    typedef typename Allocator::const_pointer   const_iterator;
    typedef typename std::reverse_iterator<iterator>       reverse_iterator;
    typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator; 

    // Accessors

//     T& operator()( int i, int j, int k, int l )
//     {
// 	return v[ index(i,j,k,l) ];
//     }
//     const T& operator()( int i, int j, int k, int l ) const
//     {
// 	return v[ index(i,j,k,l) ];
//     }

    T& operator[]( int i ) { check(i); return v[i]; }
    const T& operator[]( int i ) const { check(i); return v[i]; }

    iterator       begin()       { return v + offset(xmin,ymin,zmin,wmin); }
    const_iterator begin() const { return v + offset(xmin,ymin,zmin,wmin); }

    iterator end()
    {
	return v + offset(xmax(),ymax(),zmax(),wmax()) + 1;
    }
    const_iterator end() const
    {
	return v + offset(xmax(),ymax(),zmax(),wmax()) + 1;
    }

    reverse_iterator       rbegin()       { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(end()); }

    reverse_iterator       rend()         { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(begin()); }

//     int get_xmin() const { return xmin; }
//     int get_xlen() const { return xlen; }
//     int get_ymin() const { return ymin; }
//     int get_ylen() const { return ylen; }
//     int get_zmin() const { return zmin; }
//     int get_zlen() const { return zlen; }
//     int get_wmin() const { return wmin; }
//     int get_wlen() const { return wlen; }

    size_type size() const { return nx() * ny() * nz() * nw(); }
    size_type max_size () const { return alloc.max_size(); }
    bool empty() const { return (this->size() == 0); }

    // For backward compatibility.
    int nx() const { return xlen; }
    int ny() const { return ylen; }
    int nz() const { return zlen; }
    int nw() const { return wlen; }

    // Constructors

//     Mat4()
// 	: xmin(0), xlen(0), ymin(0), ylen(0),
// 	  zmin(0), zlen(0), wmin(0), wlen(0),
// 	  may_free_space(false), v(0)
//     {}

    Mat4( int xmax_, int ymax_, int zmax_, int wmax_, const T& t = T() )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  zmin(0), zlen(zmax_),
	  wmin(0), wlen(wmax_),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat4( T *vv, int xmax_, int ymax_, int zmax_, int wmax_ )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  zmin(0), zlen(zmax_),
	  wmin(0), wlen(wmax_),
	  may_free_space(false),
          alloc(),
          v(vv)
    {/*empty*/}

//     Mat4( const Bounds& bx, const Bounds& by,
// 	  const Bounds& bz, const Bounds& bw, const T& t = T() )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  ymin( by.min() ), ylen( by.len() ),
// 	  zmin( bz.min() ), zlen( bz.len() ),
// 	  wmin( bw.min() ), wlen( bw.len() ),
// 	  may_free_space(true),
// 	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin) )
//     {
// 	std::uninitialized_fill( begin(), end(), t );
//     }

//     Mat4( T *vv,
// 	  const Bounds& bx, const Bounds& by,
// 	  const Bounds& bz, const Bounds& bw )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  ymin( by.min() ), ylen( by.len() ),
// 	  zmin( bz.min() ), zlen( bz.len() ),
// 	  wmin( bw.min() ), wlen( bw.len() ),
// 	  may_free_space(false),
// 	  v( vv - offset(xmin,ymin,zmin,wmin) )
//     {}

    Mat4( const Mat4<T>& m )
	: xmin(m.xmin), xlen(m.xlen),
	  ymin(m.ymin), ylen(m.ylen),
	  zmin(m.zmin), zlen(m.zlen),
	  wmin(m.wmin), wlen(m.wlen),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin) )
    {
	std::uninitialized_copy( m.begin(), m.end(), begin() );
    }

    // Destructor

    ~Mat4()
    {
	detach();
    }

    void swap ( Mat4& m )
    {
        int itemp;
        bool btemp;
        Allocator atemp;
        T* ptemp;

        itemp = xmin;
        xmin = m.xmin;
        m.xmin = itemp;

        itemp = xlen;
        xlen = m.xlen;
        m.xlen = itemp;

        itemp = ymin;
        ymin = m.ymin;
        m.ymin = itemp;

        itemp = ylen;
        ylen = m.ylen;
        m.ylen = itemp;

        itemp = zmin;
        zmin = m.zmin;
        m.zmin = itemp;

        itemp = zlen;
        zlen = m.zlen;
        m.zlen = itemp;

        itemp = wmin;
        wmin = m.wmin;
        m.wmin = itemp;

        itemp = wlen;
        wlen = m.wlen;
        m.wlen = itemp;

        btemp = may_free_space;
        may_free_space = m.may_free_space;
        m.may_free_space = btemp;

        atemp = alloc;
        alloc = m.alloc;
        m.alloc = atemp;

        ptemp = v;
        v = m.v;
        m.v = ptemp;
    }

    // Boolean operators

    bool operator==( const Mat4& m ) const
    {
        if (this == &m)
            return true;

        if ( m.size() != this->size() )
            return false;

        const_iterator miter = m.begin();
        for (const_iterator iter = this->begin(); iter != this->end(); ++iter)
        {
            if (*iter != *miter) return false;
            ++miter;
        }

        return true;
    }

    bool operator!=( const Mat4& m ) const
    {
        return !(*this == m);
    }

    bool operator<( const Mat4& m ) const
    {
        return std::lexicographical_compare(this->begin(), this->end(), 
					    m.begin(), m.end());
    }

    bool operator>( const Mat4& m ) const
    {
        return (m < *this);
    }

    bool operator<=( const Mat4& m ) const
    {
        return !(m < *this);
    }

    bool operator>=( const Mat4& m ) const
    {
        return !(*this < m);
    }

    // Assignment operators

//     Mat4& operator=( const T& t )
//     {
// 	std::fill( begin(), end(), t );
// 	return *this;
//     }

    Mat4& operator=( const Mat4& m )
    {
	if (this == &m) return *this;

	if ( m.xmin != xmin || m.xlen != xlen ||
	     m.ymin != ymin || m.ylen != ylen ||
	     m.zmin != zmin || m.zlen != zlen ||
	     m.wmin != wmin || m.wlen != wlen ) {
	    detach();
	    xmin = m.xmin;
	    xlen = m.xlen;
	    ymin = m.ymin;
	    ylen = m.ylen;
	    zmin = m.zmin;
	    zlen = m.zlen;
	    wmin = m.wmin;
	    wlen = m.wlen;
	    v = alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin);
	    std::uninitialized_copy( m.begin(), m.end(), begin() );
	    may_free_space = true;
	}
	else {
	    if (v)
		std::copy( m.begin(), m.end(), begin() );
	}

	return *this;
    }

    // Mathematical support

//     template<class X> Mat4& operator+=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ += x;
// 	return *this;
//     }
//     Mat4& operator+=( const Mat4<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ += *j++;
// 	return *this;
//     }

//     template<class X> Mat4& operator-=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ -= x;
// 	return *this;
//     }
//     Mat4& operator-=( const Mat4<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ -= *j++;
// 	return *this;
//     }

//     template<class X> Mat4& operator*=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ *= x;
// 	return *this;
//     }
//     Mat4& operator*=( const Mat4<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ *= *j++;
// 	return *this;
//     }

//     template<class X> Mat4& operator/=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ /= x;
// 	return *this;
//     }
//     Mat4& operator/=( const Mat4<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ /= *j++;
// 	return *this;
//     }

    // Utility support

//     void assert_conformality( const Mat4<T>& m ) const
//     {
// 	Assert( xmin == m.xmin );
// 	Assert( xlen == m.xlen );
// 	Assert( ymin == m.ymin );
// 	Assert( ylen == m.ylen );
// 	Assert( zmin == m.zmin );
// 	Assert( zlen == m.zlen );
// 	Assert( wmin == m.wmin );
// 	Assert( wlen == m.wlen );
//     }

//     void redim( int nxmax, int nymax,
// 		int nzmax, int nwmax, const T& t = T() )
//     {
// 	// This one only works right if xmin == 0 and ymin == 0 and zmin == 0 and
// 	// wmin == 0.
// 	Assert( xmin == 0 );
// 	Assert( ymin == 0 );
// 	Assert( zmin == 0 );
// 	Assert( wmin == 0 );
// 	if (v && !may_free_space) {
// 	    // User thinks he wants to expand the aliased region.
// 	    xlen = nxmax;
// 	    ylen = nymax;
// 	    zlen = nzmax;
// 	    wlen = nwmax;
// 	    return;
// 	}
// 	detach();
// 	xlen = nxmax;
// 	ylen = nymax;
// 	zlen = nzmax;
// 	wlen = nwmax;
// 	v = alloc.allocate( size() ) - offset(0,0,0,0);
// 	std::uninitialized_fill( begin(), end(), t );
// 	may_free_space = true;
//     }

//     void redim( const Bounds& bx, const Bounds& by,
// 		const Bounds& bz, const Bounds& bw, const T& t = T() )
//     {
// 	if (v && !may_free_space) {
// 	    // Respecify the aliased region.
// 	    v += offset(xmin,ymin,zmin,wmin);
// 	    xmin = bx.min(); xlen = bx.len();
// 	    ymin = by.min(); ylen = by.len();
// 	    zmin = bz.min(); zlen = bz.len();
// 	    wmin = bw.min(); wlen = bw.len();
// 	    v -= offset(xmin,ymin,zmin,wmin);
// 	    return;
// 	}
// 	detach();
// 	xmin = bx.min(); xlen = bx.len();
// 	ymin = by.min(); ylen = by.len();
// 	zmin = bz.min(); zlen = bz.len();
// 	wmin = bw.min(); wlen = bw.len();
// 	v = alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin);
// 	std::uninitialized_fill( begin(), end(), t );
// 	may_free_space = true;
//     }

    // WARNING: This doesn't make a lot of sense anymore.

    // Check to see if this Mat4<T> is of size x by y by z by w.

//     bool conformal( int x, int y, int z, int w ) const
//     {
// 	return x == xlen && y == ylen && z == zlen && w == wlen;
//     }

//     // Obtain dimensions of this Mat4<T>.

//     void elements( int& nx, int& ny, int& nz, int& nw ) const
//     {
// 	nx = xlen; ny = ylen; nz = zlen; nw = wlen;
//     }
};

//===========================================================================//
/*!
 * \class Mat5
 * \brief A 5-d container.
 *
 * Mat5 is a container which supports four indices.  It is otherwise similar 
 * to Mat1, Mat2, Mat3 and Mat4.
 */
/*!
 * \example ds++/test/tstMat5RA.cc
 * 
 * Test of Mat5.
 */
//===========================================================================//

template< typename T,
	  typename Allocator = typename std::allocator<T> >
class Mat5
{
  private:
    int xmin, xlen, ymin, ylen, zmin, zlen, wmin, wlen, umin, ulen;
    bool may_free_space;

    int xmax() const { return xmin + xlen - 1; }
    int ymax() const { return ymin + ylen - 1; }
    int zmax() const { return zmin + zlen - 1; }
    int wmax() const { return wmin + wlen - 1; }
    int umax() const { return umin + ulen - 1; }

    // index() is used for indexing, so is checked.  offset() is used for
    // adjusting the memory pointer, which is logically distinct from
    // indexing, so is not checked.  Note in particular, that a zero size
    // matrix will have no valid index, but we still need to be able to
    // compute the offsets.

    // Compute the offset into the data array, of the i,j,k,l,m th element.
//     int index( int i, int j, int k, int l, int m ) const
//     {
// 	Assert( i >= xmin );
// 	Assert( i <  xmin + xlen );
// 	Assert( j >= ymin );
// 	Assert( j <  ymin + ylen );
// 	Assert( k >= zmin );
// 	Assert( k <  zmin + zlen );
// 	Assert( l >= wmin );
// 	Assert( l <  wmin + wlen );
// 	Assert( m >= umin );
// 	Assert( m <  umin + ulen );

// 	return offset(i,j,k,l,m);
//     }
    int offset( int i, int j, int k, int l, int m ) const
    {
	return xlen*(ylen*(zlen*(wlen*m+l)+k)+j)+i;
    }

    // Make sure a bare integer index is within the appropriate range.
    void check( int Remember(i) ) const
    {
	Assert( i >= offset( xmin,   ymin,   zmin,   wmin,   umin   ) );
	Assert( i <= offset( xmax(), ymax(), zmax(), wmax(), umax() ) );
    }

    void detach()
    {
	if (may_free_space) {
	    rtt_dsxx::Destroy( begin(), end() );
	    alloc.deallocate( v + offset(xmin,ymin,zmin,wmin,umin), size() );
	}
    }

  protected:
    Allocator alloc;
    T *v;

  public:
    typedef       T  value_type;
    typedef       T& reference;
    typedef const T& const_reference;
    typedef       T* pointer;
    typedef const T* const_pointer;
    typedef typename Allocator::difference_type difference_type;
    typedef typename Allocator::size_type       size_type;
    typedef typename Allocator::pointer         iterator;
    typedef typename Allocator::const_pointer   const_iterator;
    typedef typename std::reverse_iterator<iterator>       reverse_iterator;
    typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;

    // Accessors

//     T& operator()( int i, int j, int k, int l, int m )
//     {
// 	return v[ index(i,j,k,l,m) ];
//     }
//     const T& operator()( int i, int j, int k, int l, int m ) const
//     {
// 	return v[ index(i,j,k,l,m) ];
//     }

    T& operator[]( int i ) { check(i); return v[i]; }
    const T& operator[]( int i ) const { check(i); return v[i]; }

    iterator       begin()       { return v + offset(xmin,ymin,zmin,wmin,umin); }
    const_iterator begin() const { return v + offset(xmin,ymin,zmin,wmin,umin); }

    iterator end()
    {
	return v + offset(xmax(),ymax(),zmax(),wmax(),umax()) + 1;
    }
    const_iterator end() const
    {
	return v + offset(xmax(),ymax(),zmax(),wmax(),umax()) + 1;
    }

    reverse_iterator       rbegin()       { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(end()); }

    reverse_iterator       rend()         { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(begin()); }

//     int get_xmin() const { return xmin; }
//     int get_xlen() const { return xlen; }
//     int get_ymin() const { return ymin; }
//     int get_ylen() const { return ylen; }
//     int get_zmin() const { return zmin; }
//     int get_zlen() const { return zlen; }
//     int get_wmin() const { return wmin; }
//     int get_wlen() const { return wlen; }
//     int get_umin() const { return umin; }
//     int get_ulen() const { return ulen; }

    size_type size() const { return nx() * ny() * nz() * nw() * nu(); }
    size_type max_size () const { return alloc.max_size(); }
    bool empty() const { return (this->size() == 0); }

    // For backward compatibility.
    int nx() const { return xlen; }
    int ny() const { return ylen; }
    int nz() const { return zlen; }
    int nw() const { return wlen; }
    int nu() const { return ulen; }

    // Constructors

//     Mat5()
// 	: xmin(0), xlen(0), ymin(0), ylen(0),
// 	  zmin(0), zlen(0), wmin(0), wlen(0), umin(0), ulen(0),
// 	  may_free_space(false), v(0)
//     {}

    Mat5( int xmax_, int ymax_, int zmax_, int wmax_, int umax_,
	  const T& t = T() )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  zmin(0), zlen(zmax_),
	  wmin(0), wlen(wmax_),
	  umin(0), ulen(umax_),
	  may_free_space(true),
	  alloc(),
          v( alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin,umin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat5( T *vv, int xmax_, int ymax_, int zmax_, int wmax_, int umax_ )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  zmin(0), zlen(zmax_),
	  wmin(0), wlen(wmax_),
	  umin(0), ulen(umax_),
	  may_free_space(false),
          alloc(),
          v(vv)
    {/*empty*/}

//     Mat5( const Bounds& bx, const Bounds& by,
// 	  const Bounds& bz, const Bounds& bw,
// 	  const Bounds& bu, const T& t = T() )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  ymin( by.min() ), ylen( by.len() ),
// 	  zmin( bz.min() ), zlen( bz.len() ),
// 	  wmin( bw.min() ), wlen( bw.len() ),
// 	  umin( bu.min() ), ulen( bu.len() ),
// 	  may_free_space(true),
// 	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin,umin) )
//     {
// 	std::uninitialized_fill( begin(), end(), t );
//     }

//     Mat5( T *vv,
// 	  const Bounds& bx, const Bounds& by,
// 	  const Bounds& bz, const Bounds& bw, const Bounds& bu )
// 	: xmin( bx.min() ), xlen( bx.len() ),
// 	  ymin( by.min() ), ylen( by.len() ),
// 	  zmin( bz.min() ), zlen( bz.len() ),
// 	  wmin( bw.min() ), wlen( bw.len() ),
// 	  umin( bu.min() ), ulen( bu.len() ),
// 	  may_free_space(false),
// 	  v( vv - offset(xmin,ymin,zmin,wmin,umin) )
//     {}

    Mat5( const Mat5<T>& m )
	: xmin(m.xmin), xlen(m.xlen),
	  ymin(m.ymin), ylen(m.ylen),
	  zmin(m.zmin), zlen(m.zlen),
	  wmin(m.wmin), wlen(m.wlen),
	  umin(m.umin), ulen(m.ulen),
	  may_free_space(true),
          alloc(),
	  v( alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin,umin) )
    {
	std::uninitialized_copy( m.begin(), m.end(), begin() );
    }

    // Destructor

    ~Mat5()
    {
	detach();
    }

    void swap ( Mat5& m )
    {
        int itemp;
        bool btemp;
        Allocator atemp;
        T* ptemp;

        itemp = xmin;
        xmin = m.xmin;
        m.xmin = itemp;

        itemp = xlen;
        xlen = m.xlen;
        m.xlen = itemp;

        itemp = ymin;
        ymin = m.ymin;
        m.ymin = itemp;

        itemp = ylen;
        ylen = m.ylen;
        m.ylen = itemp;

        itemp = zmin;
        zmin = m.zmin;
        m.zmin = itemp;

        itemp = zlen;
        zlen = m.zlen;
        m.zlen = itemp;

        itemp = wmin;
        wmin = m.wmin;
        m.wmin = itemp;

        itemp = wlen;
        wlen = m.wlen;
        m.wlen = itemp;

        itemp = umin;
        umin = m.umin;
        m.umin = itemp;

        itemp = ulen;
        ulen = m.ulen;
        m.ulen = itemp;

        btemp = may_free_space;
        may_free_space = m.may_free_space;
        m.may_free_space = btemp;

        atemp = alloc;
        alloc = m.alloc;
        m.alloc = atemp;

        ptemp = v;
        v = m.v;
        m.v = ptemp;
    }

    // Boolean operators

    bool operator==( const Mat5& m ) const
    {
        if (this == &m)
            return true;

        if ( m.size() != this->size() )
            return false;

        const_iterator miter = m.begin();
        for (const_iterator iter = this->begin(); iter != this->end(); ++iter)
        {
            if (*iter != *miter) return false;
            ++miter;
        }

        return true;
    }

    bool operator!=( const Mat5& m ) const
    {
        return !(*this == m);
    }

    bool operator<( const Mat5& m ) const
    {
        return std::lexicographical_compare(this->begin(), this->end(), 
					    m.begin(), m.end());
    }

    bool operator>( const Mat5& m ) const
    {
        return (m < *this);
    }

    bool operator<=( const Mat5& m ) const
    {
        return !(m < *this);
    }

    bool operator>=( const Mat5& m ) const
    {
        return !(*this < m);
    }

    // Assignment operators

//     Mat5& operator=( const T& t )
//     {
// 	std::fill( begin(), end(), t );
// 	return *this;
//     }

    Mat5& operator=( const Mat5& m )
    {
	if (this == &m) return *this;

	if ( m.xmin != xmin || m.xlen != xlen ||
	     m.ymin != ymin || m.ylen != ylen ||
	     m.zmin != zmin || m.zlen != zlen ||
	     m.wmin != wmin || m.wlen != wlen ||
	     m.umin != umin || m.ulen != ulen ) {
	    detach();
	    xmin = m.xmin;
	    xlen = m.xlen;
	    ymin = m.ymin;
	    ylen = m.ylen;
	    zmin = m.zmin;
	    zlen = m.zlen;
	    wmin = m.wmin;
	    wlen = m.wlen;
	    umin = m.umin;
	    ulen = m.ulen;
	    v = alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin,umin);
	    std::uninitialized_copy( m.begin(), m.end(), begin() );
	    may_free_space = true;
	}
	else {
	    if (v)
		std::copy( m.begin(), m.end(), begin() );
	}

	return *this;
    }

    // Mathematical support

//     template<class X> Mat5& operator+=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ += x;
// 	return *this;
//     }
//     Mat5& operator+=( const Mat5<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ += *j++;
// 	return *this;
//     }

//     template<class X> Mat5& operator-=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ -= x;
// 	return *this;
//     }
//     Mat5& operator-=( const Mat5<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ -= *j++;
// 	return *this;
//     }

//     template<class X> Mat5& operator*=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ *= x;
// 	return *this;
//     }
//     Mat5& operator*=( const Mat5<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ *= *j++;
// 	return *this;
//     }

//     template<class X> Mat5& operator/=( const X& x )
//     {
// 	for( iterator i = begin(); i != end(); )
// 	    *i++ /= x;
// 	return *this;
//     }
//     Mat5& operator/=( const Mat5<T>& m )
//     {
// 	assert_conformality( m );
// 	iterator i = begin();
// 	const_iterator j = m.begin();
// 	while( i != end() ) *i++ /= *j++;
// 	return *this;
//     }

    // Utility support

//     void assert_conformality( const Mat5<T>& m ) const
//     {
// 	Assert( xmin == m.xmin );
// 	Assert( xlen == m.xlen );
// 	Assert( ymin == m.ymin );
// 	Assert( ylen == m.ylen );
// 	Assert( zmin == m.zmin );
// 	Assert( zlen == m.zlen );
// 	Assert( wmin == m.wmin );
// 	Assert( wlen == m.wlen );
// 	Assert( umin == m.umin );
// 	Assert( ulen == m.ulen );
//     }

//     void redim( int nxmax, int nymax, int nzmax,
// 		int nwmax, int numax, const T& t = T() )
//     {
// 	// This one only works right if xmin == 0 and ymin == 0 and zmin == 0
// 	// and wmin == 0 and umin == 0.
// 	Assert( xmin == 0 );
// 	Assert( ymin == 0 );
// 	Assert( zmin == 0 );
// 	Assert( wmin == 0 );
// 	Assert( umin == 0 );
// 	if (v && !may_free_space) {
// 	    // User thinks he wants to expand the aliased region.
// 	    xlen = nxmax;
// 	    ylen = nymax;
// 	    zlen = nzmax;
// 	    wlen = nwmax;
// 	    ulen = numax;
// 	    return;
// 	}
// 	detach();
// 	xlen = nxmax;
// 	ylen = nymax;
// 	zlen = nzmax;
// 	wlen = nwmax;
// 	ulen = numax;
// 	v = alloc.allocate( size() ) - offset(0,0,0,0,0);
// 	std::uninitialized_fill( begin(), end(), t );
// 	may_free_space = true;
//     }

//     void redim( const Bounds& bx, const Bounds& by, const Bounds& bz, 
// 		const Bounds& bw, const Bounds& bu, const T& t = T() )
//     {
// 	if (v && !may_free_space) {
// 	    // Respecify the aliased region.
// 	    v += offset(xmin,ymin,zmin,wmin,umin);
// 	    xmin = bx.min(); xlen = bx.len();
// 	    ymin = by.min(); ylen = by.len();
// 	    zmin = bz.min(); zlen = bz.len();
// 	    wmin = bw.min(); wlen = bw.len();
// 	    umin = bu.min(); ulen = bu.len();
// 	    v -= offset(xmin,ymin,zmin,wmin,umin);
// 	    return;
// 	}
// 	detach();
// 	xmin = bx.min(); xlen = bx.len();
// 	ymin = by.min(); ylen = by.len();
// 	zmin = bz.min(); zlen = bz.len();
// 	wmin = bw.min(); wlen = bw.len();
// 	umin = bu.min(); ulen = bu.len();
// 	v = alloc.allocate( size() ) - offset(xmin,ymin,zmin,wmin,umin);
// 	std::uninitialized_fill( begin(), end(), t );
// 	may_free_space = true;
//     }

    // WARNING: This doesn't make a lot of sense anymore.

    // Check to see if this Mat5<T> is of size x by y by z by w by u.

//     bool conformal( int x, int y, int z, int w, int u ) const
//     {
// 	return x == xlen && y == ylen && z == zlen
// 	    && w == wlen && u == ulen;
//     }

    // Obtain dimensions of this Mat5<T>.

//     void elements( int& nx, int& ny, int& nz, int& nw, int& nu ) const
//     {
// 	nx = xlen; ny = ylen; nz = zlen; nw = wlen; nu = ulen;
//     }
};

} // end of rtt_dsxx

#endif                          // rtt_ds_Mat_hh

//---------------------------------------------------------------------------//
//                              end of ds++/Mat.hh
//---------------------------------------------------------------------------//
