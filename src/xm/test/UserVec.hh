//----------------------------------*-C++-*----------------------------------//
/*! \file    UserVec.hh
 *  \author  Geoffrey Furnish
 *  \date    Thu Sep 12 10:51:22 1996
 *  \brief
 *  \note    Copyright (C) 1996-2010 Los Alamos National Security, LLC.
 *           All rights reserved.
 *  \version $Id$
 */
//---------------------------------------------------------------------------//

#ifndef __xm_UserVec_hh__
#define __xm_UserVec_hh__

#include "../xm.hh"
#include "ds++/Assert.hh"

//===========================================================================//
/*! \class UserVec - Prototypical "user" array class
 *
 * The purpose of this class is to serve as a prototypical example of a user
 * defined indexable class, which is to be integrated into the expression math
 * suite.  The goal is very simple: to see how a user's indexable class may be
 * included into the expression math suite, without modifying the expression
 * math suite.  That is, we want to find out what a client can do to include
 * his own indexable classes in expression math statements, by making only
 * reasonable modifications to his own class, without touching the xm suite
 * itself.  Here goes...
 */
//===========================================================================//
template<class T>
class UserVec : public xm::Indexable<T,UserVec<T> >
{

    T *v;
    int sz;

    UserVec( const UserVec& );

  public:
    explicit UserVec( int n )
	: xm::Indexable<T,UserVec<T> >(),
          v(NULL),
	  sz(n) 
    {
	v = new T[sz];
    }

    ~UserVec() { delete[] v; }

    UserVec<T>& operator=( T t )
    {
	for( int i=0; i < sz; i++ ) v[i] = t;
	return *this;
    }

    UserVec& operator=( const UserVec& x)
    {
	for( int i=0; i < sz; i++ ) v[i] = x[i];
	return *this;
    }

    int size() const { return sz; }

    T& operator[]( int n ) {
	Assert( n >= 0 && n < sz );
	return v[n];
    }
    T operator[]( int n ) const {
	Assert( n >= 0 && n < sz );
	return v[n];
    }

    // This is the routine which evaluates the expression.
    template<class X>
    UserVec<T>& operator=( const xm::Xpr< T, X, UserVec<T> >& x )
    {
	// we need to scope the Base class because some compilers (CXX 6.5.0)
	// suck
 	return xm::Indexable<T,UserVec<T> >::assign_from( x );
    }
};

//===========================================================================//
// class FooBar - A test class.
//===========================================================================//

template<class T>
class FooBar : public xm::Indexable<T,FooBar<T> >
{
    T *v;
    int sz;

    FooBar( const FooBar& );
    FooBar& operator=( const FooBar& );

  public:
    FooBar( size_t const n )
	: xm::Indexable<T,FooBar<T> >(),
          v(NULL),
	  sz(n) 
    {
	v = new T[sz];
    }

    ~FooBar() { delete[] v; }

    FooBar<T>& operator=( T t )
    {
	for( int i=0; i < sz; i++ ) v[i] = t;
	return *this;
    }

    int size() const { return sz; }

    T& operator[]( int n ) {
	Assert( n >= 0 && n < sz );
	return v[n];
    }
    T operator[]( int n ) const {
	Assert( n >= 0 && n < sz );
	return v[n];
    }

// This is the routine which evaluates the expression.

    template<class X>
    FooBar<T>& operator=( const xm::Xpr< T, X, FooBar<T> >& x )
    {
	// we need to scope the Base class because some compilers (CXX 6.5.0)
	// suck
 	return xm::Indexable<T,FooBar<T> >::assign_from( x );
    }
};

#endif // __xm_UserVec_hh__

//---------------------------------------------------------------------------//
// end of UserVec.hh
//---------------------------------------------------------------------------//
