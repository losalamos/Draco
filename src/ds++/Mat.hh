//----------------------------------*-C++-*----------------------------------//
// Mat.hh
// Geoffrey Furnish
// Fri Jan 24 15:48:31 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __ds_Mat_hh__
#define __ds_Mat_hh__

// The matrix family of classes presented in this file represent a total
// rewrite of the matrix family used in GTS up until this time.  The C++
// industry has matured dramatically since the initial days of GTS.  There
// is now a pretty decent, and very thoroughly documented, standard library
// for C++, which largely subsums the original goals of DS++.  Moreover, the
// compiler industry has finally begun to produce compilers that are
// acceptable in terms of their language conformance.  This necesitates some
// "repositioning" of DS++.  IMO, people don't really need a map class
// anymore, since they can use the one in the STL.  But there are things the
// STL doesn't really cover, and that is where DS++ is going to have to be
// positioned in order to be useful.  The first effort at modernizing DS++
// will thus concentrate on moving the Matrix family out of GTS proper and
// into DS++.  This will involve generalizing it so that it works more like
// an STL container.  The reconstituted family will be just a container
// family.  The algorithms which were member functions in the old GTS matrix
// family will be split out, just as in the STL.  What we will retain,
// however, and improve in the process, is the ability to provide a
// sophisticated container facility with good memory management.  Also, we
// will fix the instantiation annoyances which have haunted the older GTS
// matrix family.

#include <algorithm>

#include "Assert.hh"
#include "Allocators.hh"
#include "Bounds.hh"

#ifdef TEST_MAT
#define private public
#define protected public
#endif

//===========================================================================//
// class Mat1 - A 1-d container.

// This class is intended to be a replacement for the STL vector<T> class.
// The reasons for wanting a replacement for vector<T> are primarily
// threefold: 1) Integration with the DS++ assertion model.  2) Ability to
// provide more sophisticated allocators, and 3) provide a 1-d analog to Mat2
// and Mat3.
//===========================================================================//
#if 0
template< class T, class Allocator = typename alloc_traits<T>::Default_Allocator >
class Mat1 {
    int xmin, xmax;
    bool may_free_space;

    void check( int i ) const { Assert( i >= xmin && i < xmax ); }

  protected:
    Allocator alloc;
    T *v;

  public:
    typedef typename Allocator::iterator iterator;
    typedef typename Allocator::const_iterator const_iterator;

// Constructors

    Mat1() : xmin(0), xmax(0), may_free_space(false), v(0) {}

    Mat1( int _xmax, const T& t = T() )
	: xmin(0), xmax(_xmax), may_free_space(true),
	  v(alloc.fetch(xmax))
    {
	std::uninitialized_fill_n( v+xmin, nx(), t );
    }
    
    Mat1( const Mat1<T>& m )
	: xmin(m.xmin), xmax(m.xmax),
	  may_free_space(true), v(alloc.fetch(xmax-xmin) - xmin)
    {
	std::copy( m.begin(), m.end(), v+xmin );
    }

    Mat1( T *vv, int x )
	: xmin(0), xmax(x), may_free_space(false), v(vv) {}

    Mat1( T *begin, T *end )
	: xmin(0), xmax(end - begin),
	  may_free_space(false), v(begin)
    {}

    Mat1( const Bounds& b, const T& t = T() )
	: xmin( b.min() ), xmax( b.max() + 1 ),
	  may_free_space(true),
	  v( alloc.fetch(xmax-xmin) - xmin )
    {
	std::uninitialized_fill_n( v+xmin, nx(), t );
    }

    Mat1( T *vv, const Bounds& b )
	: xmin( b.min() ), xmax( b.max() + 1 ),
	  may_free_space(false),
	  v(vv)
    {}

// Destructor

    ~Mat1()
    {
	if (may_free_space) {
	    std::destroy( v+xmin, v+xmax );
	    alloc.release( v+xmin, xmax-xmin );
	}
    }

// Assignment operators

    Mat1& operator=( const T& t )
    {
	for( int i=xmin; i < xmax; i++ )
	    v[i] = t;

	return *this;
    }

    Mat1& operator=( const Mat1& m )
    {
	if (this == &m) return *this;

	if (m.xmin != xmin || m.xmax != xmax) {
	    if (may_free_space) alloc.release( v+xmin, xmax-xmin );
	    xmin = m.xmin;
	    xmax = m.xmax;
	    v = alloc.fetch( xmax ) - xmin;
	}

	std::copy( m.begin(), m.end(), v+xmin );

	return *this;
    }

// Accessors

    T& operator()( int i ) { check(i); return v[i]; }
    T& operator[]( int i ) { check(i); return v[i]; }

    const T& operator()( int i ) const { check(i); return v[i]; }
    const T& operator[]( int i ) const { check(i); return v[i]; }

    iterator begin() { return v+xmin; }
    const_iterator begin() const { return v+xmin; }

    iterator end() { return v+xmax; }
    const_iterator end() const { return v+xmax; }

    int nx()   const { return xmax - xmin; }
    int size() const { return xmax - xmin; }
    int base() const { return xmin; }

// Mathematical support

    template<class X> Mat1& operator+=( const X& x )
    {
	for( int i=xmin; i < xmax; i++ )
	    v[i] += x;
	return *this;
    }
    Mat1& operator+=( const Mat1<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	for( int i=xmin; i < xmax; i++ )
	    v[i] += m.v[i];
	return *this;
    }

    template<class X> Mat1& operator-=( const X& x )
    {
	for( int i=xmin; i < xmax; i++ )
	    v[i] -= x;
	return *this;
    }
    Mat1& operator-=( const Mat1<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	for( int i=xmin; i < xmax; i++ )
	    v[i] -= m.v[i];
	return *this;
    }

    template<class X> Mat1& operator*=( const X& x )
    {
	for( int i=xmin; i < xmax; i++ )
	    v[i] *= x;
	return *this;
    }
    Mat1& operator*=( const Mat1<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	for( int i=xmin; i < xmax; i++ )
	    v[i] *= m.v[i];
	return *this;
    }

    template<class X> Mat1& operator/=( const X& x )
    {
	for( int i=xmin; i < xmax; i++ )
	    v[i] /= x;
	return *this;
    }
    Mat1& operator/=( const Mat1<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	for( int i=xmin; i < xmax; i++ )
	    v[i] /= m.v[i];
	return *this;
    }

// Utility support

    void redim( int nxmax, const T& t = T() )
    {
	if (v && !may_free_space) {
	// User thinks he wants to expand the aliased region.
	    xmax = nxmax;
	    return;
	}
	if (v) {
	    std::destroy( v+xmin, v+xmax );
	    alloc.release( v+xmin, xmax-xmin );
	}
	xmax = nxmax;
	v = alloc.fetch(xmax-xmin) - xmin;
	may_free_space = true;
	std::uninitialized_fill_n( v, nx(), t );
    }

// WARNING: This doesn't make a lot of sense anymore.
    bool conformal( int n ) const { return n == xmax; }
};
#endif

template< class T,
    class Allocator = typename alloc_traits<T>::Default_Allocator >
class Mat1 {
  private:
    int xmin, xlen;
    bool may_free_space;

    int xmax() const { return xmin + xlen - 1; }

    int index( int i ) const
    {
	Assert( i >= xmin );
	Assert( i < xmin + xlen );
	return i;
    }

    void detach()
    {
	if (may_free_space) {
	    std::destroy( begin(), end() );
	    alloc.release( v + index(xmin), size() );
	}
    }

  protected:
    Allocator alloc;
    T *v;

    T&       operator[]( int i )       { return v[ index(i) ]; }
    const T& operator[]( int i ) const { return v[ index(i) ]; }

  public:
    typedef typename Allocator::iterator iterator;
    typedef typename Allocator::const_iterator const_iterator;

// Accessors

    T&       operator()( int i )       { return v[ index(i) ]; }
    const T& operator()( int i ) const { return v[ index(i) ]; }

    iterator       begin()       { return v + index(xmin); }
    const_iterator begin() const { return v + index(xmin); }

    iterator       end()         { return v + index(xmax) + 1; }
    const_iterator end() const   { return v + index(xmax) + 1; }

    int nx() const { return xlen; }
    int size() const { return nx(); }

// Constructors

    Mat1()
	: xmin(0), xlen(0),
	  may_free_space(false), v(0)
    {}

    Mat1( int xmax_, const T& t = T() )
	: xmin(0), xlen(xmax_),
	  may_free_space(true),
	  v( alloc.fetch( size() ) - index(xmin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat1( T *vv, int xmax_ )
	: xmin(0), xlen(xmax_),
	  may_free_space(false), v(vv)
    {}

    Mat1( const Bounds& bx, const T& t = T() )
	: xmin( bx.min() ), xlen( bx.len() ),
	  may_free_space(true),
	  v( alloc.fetch( size() ) - index(xmin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat1( T *vv, const Bounds& bx )
	: xmin( bx.min() ), xlen( bx.len() ),
	  may_free_space(true),
	  v( vv - index(xmin) )
    {}

    Mat1( const Mat1<T>& m )
	: xmin(m.xmin), xlen(m.xlen),
	  may_free_space(true),
	  v( alloc.fetch( size() ) - index(xmin) )
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
	std::fill( begin(), end(), t )
	return *this;
    }

    Mat1& operator=( const Mat1& m )
    {
	if (this == &m) return *this;

	if ( m.xmin != xmin || m.xlen != xlen ) {
	    detach();
	    xmin = m.xmin;
	    xlen = m.xlen;
	    v = alloc.fetch( size() ) - index(xmin);
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

    template<class X> Mat1& operator+=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ += x;
	return *this;
    }
    Mat1& operator+=( const Mat1<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ += *j++;
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
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ -= *j++;
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
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ *= *j++;
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
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ /= *j++;
	return *this;
    }

// Utility support

    void redim( int nxmax, const T& t = T() )
    {
    // This one only works right if xmin == 0.
	Assert( xmin == 0 );
	if (v && !may_free_space) {
	// User thinks he wants to expand the aliased region.
	    xmax = nxmax;
	    return;
	}
	detach();
	xlen = nxmax;
	v = alloc.fetch( size() ) - index(xmin);
	std::uninitialized_fill( begin(), end(), t );
	may_free_space = true;
    }

    void redim( const Bounds& bx, const T& t = T() )
    {
	if (v && !may_free_space) {
	// Respecify the aliased region.
	    v += index(xmin);
	    xmin = bx.min(); xlen = bx.len();
	    v -= index(xmin);
	    return;
	}
	detach();
	xmin = bx.min(); xlen = bx.len();
	v = alloc.fetch( size() ) - index(xmin);
	std::uninitialized_fill( begin(), end(), t );
	may_free_space = true;
    }

// Check to see if this Mat2<T> is of size x by y.

    bool conformal( int x ) const
    {
	return xmin == 0 && x == xlen;
    }

// Obtain dimensions of this Mat2<T>.

    void elements( int& nx_ ) const { nx_ = nx(); }
};

//===========================================================================//
// class Mat2 - A 2-d container.

//===========================================================================//

template< class T, class Allocator = Simple_Allocator<T> >
class Mat2 {
  private:
    int xmin, xlen, ymin, ylen;
    bool may_free_space;

    int xmax() const { return xmin + xlen - 1; }
    int ymax() const { return ymin + ylen - 1; }

// Compute the offset into the data array, of the i,j th element.
    int index( int i, int j ) const
    {
	Assert( i >= xmin );
	Assert( i < xmin + xlen );
	Assert( j >= ymin );
	Assert( j < ymin + ylen );

	return xlen * j + i;
    }

// Make sure a bare integer index is within the appropriate range.
    void check( int i )
    {
	Assert( i >= index( xmin, ymin ) );
	Assert( i <= index( xmax, ymax ) );
    }

    void detach()
    {
	if (may_free_space) {
	    std::destroy( begin(), end() );
	    alloc.release( v +index(xmin,ymin), size() );
	}
    }

  protected:
    Allocator alloc;
    T *v;

    T& operator[]( int i ) { check(i); return v[i]; }
    const T& operator[]( int i ) const { check(i); return v[i]; }

  public:
    typedef typename Allocator::iterator iterator;
    typedef typename Allocator::const_iterator const_iterator;

// Accessors

    T&       operator()( int i, int j )       { return v[ index(i,j) ]; }
    const T& operator()( int i, int j ) const { return v[ index(i,j) ]; }

    iterator       begin()       { return v + index(xmin,ymin); }
    const_iterator begin() const { return v + index(xmin,ymin); }

    iterator       end()         { return v + index(xmax(),ymax()) + 1; }
    const_iterator end() const   { return v + index(xmax(),ymax()) + 1; }

    int nx() const { return xlen; }
    int ny() const { return ylen; }
    int size() const { return nx() * ny(); }

// Constructors

    Mat2()
	: xmin(0), xlen(0), ymin(0), ylen(0),
	  may_free_space(false), v(0)
    {}

    Mat2( int xmax_, int ymax_, const T& t = T() )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  may_free_space(true),
	  v( alloc.fetch( size() ) - index(xmin,ymin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat2( T *vv, int xmax_, int ymax_ )
	: xmin(0), xlen(xmax_),
	  ymin(0), ylen(ymax_),
	  may_free_space(false), v(vv)
    {}

    Mat2( const Bounds& bx, const Bounds& by, const T& t = T() )
	: xmin( bx.min() ), xlen( bx.len() ),
	  ymin( by.min() ), ylen( by.len() ),
	  may_free_space(true),
	  v( alloc.fetch( size() ) - index(xmin,ymin) )
    {
	std::uninitialized_fill( begin(), end(), t );
    }

    Mat2( T *vv, const Bounds& bx, const Bounds& by )
	: xmin( bx.min() ), xlen( bx.len() ),
	  ymin( by.min() ), ylen( by.len() ),
	  may_free_space(true),
	  v( vv - index(xmin,ymin) )
    {}

    Mat2( const Mat2<T>& m )
	: xmin(m.xmin), xlen(m.xlen),
	  ymin(m.ymin), ylen(m.ylen),
	  may_free_space(true),
	  v( alloc.fetch( size() ) - index(xmin,ymin) )
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
	std::fill( begin(), end(), t )
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
	    v = alloc.fetch( size() ) - index(xmin,ymin);
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

    template<class X> Mat2& operator+=( const X& x )
    {
	for( iterator i = begin(); i != end(); )
	    *i++ += x;
	return *this;
    }
    Mat2& operator+=( const Mat2<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	Assert( ymin == m.ymin );
	Assert( ymax == m.ymax );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ += *j++;
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
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
	Assert( ymin == m.ymin );
	Assert( ylen == m.ylen );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ -= *j++;
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
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
	Assert( ymin == m.ymin );
	Assert( ylen == m.ylen );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ *= *j++;
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
	Assert( xmin == m.xmin );
	Assert( xlen == m.xlen );
	Assert( ymin == m.ymin );
	Assert( ylen == m.ylen );
	for( iterator i = begin(), j = m.begin(); i != end(); )
	    *i++ /= *j++;
	return *this;
    }

// Utility support

    void redim( int nxmax, int nymax, const T& t = T() )
    {
    // This one only works right if xmin == 0 and ymin == 0.
	Assert( xmin == 0 );
	Assert( ymin == 0 );
	if (v && !may_free_space) {
	// User thinks he wants to expand the aliased region.
	    xmax = nxmax;
	    ymax = nymax;
	    return;
	}
	detach();
	xlen = nxmax;
	ylen = nymax;
	v = alloc.fetch( size() ) - index(0,0);
	std::uninitialized_fill( begin(), end(), t );
	may_free_space = true;
    }

    void redim( const Bounds& bx, const Bounds& by, const T& t = T() )
    {
	if (v && !may_free_space) {
	// Respecify the aliased region.
	    v += index(xmin,ymin);
	    xmin = bx.min(); xlen = bx.len();
	    ymin = by.min(); ylen = by.len();
	    v -= index(xmin,ymin);
	    return;
	}
	detach();
	xmin = bx.min(); xlen = bx.len();
	ymin = by.min(); ylen = by.len();
	v = alloc.fetch( size() ) - index(xmin,ymin);
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
// class Mat3 - A 3-d container.

//===========================================================================//

template< class T, class Allocator = Simple_Allocator<T> >
class Mat3 {
    int xmin, xmax, ymin, ymax, zmin, zmax;
    bool may_free_space;

    void check( int i ) const { Assert( i >= 0 && i < xmax*ymax*zmax ); }

    void check( int i, int j, int k ) const
    {
	Assert( i >= xmin && i < xmax );
	Assert( j >= ymin && j < ymax );
	Assert( k >= zmin && k < zmax );
    }

  protected:
    Allocator alloc;
    T *v;

    T& operator[]( int i ) { check(i); return v[i]; }
    const T& operator[]( int i ) const { check(i); return v[i]; }

  public:
    typedef typename Allocator::iterator iterator;
    typedef typename Allocator::const_iterator const_iterator;

// Constructors

// WARNING: Not set up for offset bounds yet...

    Mat3()
	: xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0),
	  may_free_space(false), v(0)
    {}

    Mat3( int _xmax, int _ymax, int _zmax )
	: xmin(0), xmax(_xmax),
	  ymin(0), ymax(_ymax),
	  zmin(0), zmax(_zmax),
	  may_free_space(true),
	  v( alloc.fetch( xmax * ymax * zmax ) )
    {}
    
    Mat3( const Mat3<T>& m )
	: xmin(m.xmin), xmax(m.xmax),
	  ymin(m.ymin), ymax(m.ymax),
	  zmin(m.zmin), zmax(m.zmax),
	  may_free_space(true),
	  v( alloc.fetch( xmax * ymax * zmax ) )
    {
	std::copy( m.begin(), m.end(), v );
    }

    Mat3( T *vv, int xmx, int ymx, int zmx )
	: xmin(0), xmax(xmx),
	  ymin(0), ymax(ymx),
	  zmin(0), zmax(zmx),
	  may_free_space(false), v(vv)
    {}

//     Mat1( T *vv, int x ) : xmin(0), xmax(x), may_free_space(false), v(vv) {}

//     Mat1( const Bounds& b )
// 	: xmin( b.min() ), xmax( b.max() + 1 ),
// 	  may_free_space(true),
// 	  v( alloc.fetch(xmax-xmin) - xmin )
//     {}

// Destructor

    ~Mat3() { if (may_free_space) alloc.release( v, xmax * ymax * zmax ); }

// Assignment operators

    Mat3& operator=( const T& t )
    {
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] = t;

	return *this;
    }

    Mat3& operator=( const Mat3& m )
    {
	if (this == &m) return *this;

	if ( m.xmin != xmin || m.xmax != xmax ||
	     m.ymin != ymin || m.ymax != ymax ||
	     m.zmin != zmin || m.zmax != zmax ) {
	    if (may_free_space) alloc.release( v, xmax*ymax*zmax );
	    xmin = m.xmin;
	    xmax = m.xmax;
	    ymin = m.ymin;
	    ymax = m.ymax;
	    zmin = m.zmin;
	    zmax = m.zmax;
	    v = alloc.fetch( xmax*ymax*zmax );
	}

	std::copy( m.begin(), m.end(), v );

	return *this;
    }

// Accessors

    T& operator()( int i, int j, int k )
    {
	check(i,j,k);
	return v[ xmax * (k * ymax + j) + i ];
    }

    const T& operator()( int i, int j, int k ) const
    {
	check(i,j,k);
	return v[ xmax * (k * ymax + j) + i ];
    }

    iterator begin() { return v; }
    const_iterator begin() const { return v; }

    iterator end() { return v + xmax*ymax*zmax; }
    const_iterator end() const { return v + xmax*ymax*zmax; }

    int nx() const { return xmax - xmin; }
    int ny() const { return ymax - ymin; }
    int nz() const { return zmax - zmin; }

// Mathematical support

    template<class X>
    Mat3& operator+=( const X& x )
    {
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] += x;
	return *this;
    }

    Mat3& operator+=( const Mat3<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	Assert( ymin == m.ymin );
	Assert( ymax == m.ymax );
	Assert( zmin == m.zmin );
	Assert( zmax == m.zmax );
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] += m.v[i];
	return *this;
    }

    template<class X>
    Mat3& operator-=( const X& x )
    {
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] -= x;
	return *this;
    }

    Mat3& operator-=( const Mat3<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	Assert( ymin == m.ymin );
	Assert( ymax == m.ymax );
	Assert( zmin == m.zmin );
	Assert( zmax == m.zmax );
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] -= m.v[i];
	return *this;
    }

    template<class X>
    Mat3& operator*=( const X& x )
    {
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] *= x;
	return *this;
    }

    Mat3& operator*=( const Mat3<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	Assert( ymin == m.ymin );
	Assert( ymax == m.ymax );
	Assert( zmin == m.zmin );
	Assert( zmax == m.zmax );
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] *= m.v[i];
	return *this;
    }

    template<class X>
    Mat3& operator/=( const X& x )
    {
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] /= x;
	return *this;
    }

    Mat3& operator/=( const Mat3<T>& m )
    {
	Assert( xmin == m.xmin );
	Assert( xmax == m.xmax );
	Assert( ymin == m.ymin );
	Assert( ymax == m.ymax );
	Assert( zmin == m.zmin );
	Assert( zmax == m.zmax );
	for( int i=0; i < xmax*ymax*zmax; i++ )
	    v[i] /= m.v[i];
	return *this;
    }

    int operator==( const Mat3<T>& m ) const
    {
	if (xmax != m.xmax) return 0;
	if (ymax != m.ymax) return 0;
	if (zmax != m.zmax) return 0;

	for( int i=0; i < xmax*ymax*zmax; i++ )
	    if (v[i] != m.v[i]) return 0;

	return 1;
    }

    int operator!=( const Mat3<T>& m ) const
    {
	return !(*this == m);
    }

// Utility support

    void redim( int nxmax, int nymax, int nzmax, const T& t = T() )
    {
    // This one only works right if xmin == 0.
	Assert( xmin == 0 );
	Assert( ymin == 0 );
	Assert( zmin == 0 );
	if (v && !may_free_space) {
	// User thinks he wants to expand the aliased region.
	    xmax = nxmax;
	    ymax = nymax;
	    zmax = nzmax;
	    return;
	}
	if (v) {
	    std::destroy( begin(), end() );
	    alloc.release( v, xmax*ymax*zmax );
	}
	xmax = nxmax;
	ymax = nymax;
	zmax = nzmax;
	v = alloc.fetch( xmax*ymax*zmax );
	may_free_space = true;
	std::uninitialized_fill( begin(), end(), t );
    }

// This is here temporarily until Geoff figures out where he wants it to be.

    T norm() const
    {
	T nrm = T(0.);
	int i;
	for( i=0; i < xmax*ymax*zmax; i++ )
	    nrm += v[i]*v[i];

	return sqrt(nrm/i);
    }

// WARNING: This doesn't make a lot of sense anymore.

// Check to see if this Mat3<T> is of size x by y by z.

    bool conformal( int x, int y, int z ) const
    {
	return x == xmax && y == ymax && z == zmax;
    }

// Obtain dimensions of this Mat3<T>.

    void elements( int& nx, int& ny, int& nz ) const
    {
	nx = xmax; ny = ymax; nz = zmax;
    }

};

#ifdef TEST_MAT
#undef class
#undef private
#undef protected
#endif

#endif                          // __ds_Mat_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/Mat.hh
//---------------------------------------------------------------------------//
