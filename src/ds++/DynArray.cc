//----------------------------------*-C++-*----------------------------------//
// DynArray.cc
// Geoffrey Furnish
// 28 January 1994
//---------------------------------------------------------------------------//
// @> A dynamically growing array template class.
//---------------------------------------------------------------------------//

#include "DynArray.hh"
#include "Assert.hh"

//---------------------------------------------------------------------------//
// Principal constructor.  Specify initial size. Base defaults to zero, and
// growth factor defaults to 1.5, but both can be overridden if desired.
//---------------------------------------------------------------------------//

template<class T>
DynArray<T>::DynArray( int _sz, int _base, float gf )
    : sz(_sz), base(_base), growthfactor(gf)
{
    if (sz < 1) sz = 1;
    v = new T[sz];
    v -= base;			// v[base] is first allocated element.

// For lack of a better plan, we set lo and hi refs to base.

    hiref = lowref = base;
    defval = (T) 0;
}

//---------------------------------------------------------------------------//
// ctor.  Allows specification of default value.  Note that in order to use
// this method for setting the default, you have to also specify the growth
// factor and base.  This to prevent resolution trauma if T happens to be a
// float, etc.
//---------------------------------------------------------------------------//

template<class T>
DynArray<T>::DynArray( int _sz, int _base, T dv, float gf )
    : sz(_sz), base(_base), defval(dv), growthfactor(gf)
{
    if (sz < 1) sz = 1;
    v = new T[sz];
    v -= base;

    for( int i=base; i < base + sz; i++ )
      v[i] = defval;

    hiref = lowref = base;
}

//---------------------------------------------------------------------------//
// Initialize from existing DynArray<T>.
//---------------------------------------------------------------------------//

template<class T>
DynArray<T>::DynArray( const DynArray<T>& da )
{
    defval = da.defval;
    base = da.base;
    sz = da.sz;
    lowref = da.lowref;
    hiref = da.hiref;

    v = new T[sz];
    v -= base;

    for( int i=base; i < base + sz; i++ )
	v[i] = da.v[i];
}

//---------------------------------------------------------------------------//
// Assign from existing DynArray<T>.
//---------------------------------------------------------------------------//

template<class T>
DynArray<T>& DynArray<T>::operator=( const DynArray<T>& da )
{
    if ( this == &da )
	return *this;

    v += base;
    delete[] v;

    defval = da.defval;
    base = da.base;
    sz = da.sz;
    lowref = da.lowref;
    hiref = da.hiref;

    v = new T[sz];
    v -= base;

    for( int i=base; i < base + sz; i++ )
	v[i] = da.v[i];

    return *this;
}

//---------------------------------------------------------------------------//
// Automatically expands the DynArray<T> on reference.  Note that it returns
// a reference, so be careful of the famous "dangling reference" problem.
// Basically it's up to the user to avoid hosing his code.
//---------------------------------------------------------------------------//

template<class T>
T& DynArray<T>::operator[]( int n )
{
    int i;

    if ( n < base ) {		// underflow
	int nsz = (int) ( (base + sz - n) * growthfactor);
	int nbase = base + sz - nsz; // Find the new base.

	if ( !(nbase < n) ) {
	// User picked too small of a growth factor...
	    nbase = n-1;
	    nsz = base + sz - nbase;
	}

	T *nv = new T[nsz];	// Allocate the new array.
	nv -= nbase;

	// initialize the new array.

	for( i=base; i < base + sz; i++ )
	    nv[i] = v[i];
	for( i = nbase; i < base; i++ )
	    nv[i] = defval;

	v += base;
	delete[] v;		// Cleanup.
	sz = nsz;
	base = nbase;
	v = nv;

    } else if ( n >= base + sz ) { // overflow

	int nsz = (int) ( (n - base +1) * growthfactor);

	if ( n >= base + nsz ) {
	// User picked too small of a growthfactor.
	    nsz = n - base +1;
	}

	T *nv = new T[nsz];
	nv -= base;

	for( i=base; i < base + sz; i++ )
	  nv[i] = v[i];
	for( i = base + sz; i < base + nsz; i++ )
	  nv[i] = defval;

	v += base;
	delete[] v;
	sz = nsz;
	v = nv;
    }

// This could either be taken out altogether later, or, alternatively, could
// be just left in, leaving it up to the user to compile with nullified
// assertions if he so desires.

// Probably best to leave it in a bit longer so we can make sure there aren't
// any expansion logic bugs...

    Assert( n >= base && n < base + sz );

    if ( n < lowref ) lowref = n;
    if ( n > hiref )  hiref = n;

    return v[n];
}

//---------------------------------------------------------------------------//
// operator[] const

// Use this one for a const DynArray.  In other words, this one doesn't
// expand on reference.  It also doesn't export a reference.
//---------------------------------------------------------------------------//

template<class T>
T DynArray<T>::operator[]( int n ) const
{
    Assert( n >= base && n < base + sz );

    return v[n];
}

//---------------------------------------------------------------------------//
// Equality operator.
//---------------------------------------------------------------------------//

template<class T>
int DynArray<T>::operator==( const DynArray<T>& da ) const
{
// It doesn't really matter if base/sz are not identical, since they may be
// different as a result of different access patterns forcing different
// expansion characteristics.  The only thing that really mattes is if the
// two DynArray's hold the same data, and the only clear indication of what
// they hold is determined by lowref/hiref, which show what has been
// referenced. 

    if (lowref != da.lowref) return 0;
    if ( hiref != da.hiref ) return 0;

    for( int i=lowref; i <= hiref; i++ )
	if ( !(v[i] == da.v[i]) ) return 0;

    return 1;
}

//---------------------------------------------------------------------------//
// Inequality operator.
//---------------------------------------------------------------------------//

template<class T>
int DynArray<T>::operator!=( const DynArray<T>& da ) const
{
    return !(*this == da);
}

//---------------------------------------------------------------------------//
// Output a DynArray<T> to an ostream.  Useful for diagnostics, etc.
//---------------------------------------------------------------------------//

template<class T>
ostream& operator<<( ostream& os, const DynArray<T>& d )
{
    os << "DynArray<T>: low=" << d.low() << " high=" << d.high() << endl;
    for( int i= d.low(); i <= d.high(); i++ )
	os << d[i] << ' ';
    os << endl;

    return os;
}

//---------------------------------------------------------------------------//
//                              end of DynArray.cc
//---------------------------------------------------------------------------//
