//----------------------------------*-C++-*----------------------------------//
// DynArray.hh
// Geoffrey Furnish
// 28 January 1994
//---------------------------------------------------------------------------//
// @> A dynamically growing array template class.
//---------------------------------------------------------------------------//

#ifndef __ds_DynArray_hh__
#define __ds_DynArray_hh__

#include <iosfwd>

namespace rtt_dsxx
{

//===========================================================================//
// class DynArray<T> - A dynamically growable templated array class

// This class provides an array facility which expands dynamically on
// reference so that it is always big enough to hold what you want to put on
// it, but predimensioning is not required.  When a reference is made to an
// element which is out of bounds on either end, the array is resized,
// extending in the direction of the reference.  
//
// The default subscripting mechanism exports a reference to the internal
// datum.  This is obviously dangerous, and the user is hereby forewarned.
// Beware of dangling references!  You should assume that any reference may
// invalidate the previous result.  A const form of indexing allows retrieval
// by value without dynamic extension, if such semantics are desired for a
// particular situation.
//===========================================================================//

template<class T>
class DynArray {

    T *v;
    T defval;
    int base;
    int sz;
    float growthfactor;
    int lowref, hiref;

  public:

    DynArray( int sz_ =1, int base_ =0, float gf =1.5 );
    DynArray( int sz_, int base_, T dv, float gf );
    DynArray( const DynArray<T>& da );
    DynArray<T>& operator=( const DynArray<T>& da );

    ~DynArray() { v += base; delete[] v; }

    T Get_defval() const { return defval; }
    int Get_base() const { return base; }
    int Get_size() const { return sz; }
    float Get_growthfactor() const { return growthfactor; }
    int low()  const { return lowref; }
    int high() const { return hiref; }

    void low(  int l ) { lowref = l; }
    void high( int h ) { hiref  = h; }

// This is dangerous.  Some say you shouldn't even ever return a reference
// from a method if the memory could move underneath.  Well, this is an
// array, and we need these semantics.  Just have to leave it to the user to
// be sure he doesn't give his code a dirty fill.

    T& operator[]( int n );
    T operator[]( int n ) const;

    int operator==( const DynArray<T>& da ) const;
    int operator!=( const DynArray<T>& da ) const;
};

template<class T>
std::ostream& operator<<( std::ostream& os, const DynArray<T>& d );

#define INSTANTIATE_DynArray(a) \
template class DynArray<a>; \
template std::ostream& ostream<<( std::ostream& os, const DynArray<a>& d );

} // end of rtt_dsxx

#endif				// __ds_DynArray_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/DynArray.hh
//---------------------------------------------------------------------------//
