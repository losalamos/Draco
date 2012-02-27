//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Created on: Thu Sep 12 10:51:22 1996
// Created by: Geoffrey Furnish
// Also maintained by: Randy M. Roberts and Shawn Pautz
// 
//---------------------------------------------------------------------------//

#ifndef __xm_UserVec_hh__
#define __xm_UserVec_hh__

// #include "ds++/Assert.hh"
#define Assert(x)

#include "xm/xm.hh"

//===========================================================================//
// class UserVec - Prototypical "user" array class

// The purpose of this class is to serve as a prototypical example of a user
// defined indexable class, which is to be integrated into the expression
// math suite.  The goal is very simple:  to see how a user's indexable class
// may be included into the expression math suite, without modifying the
// expression math suite.  That is, we want to find out what a client can do
// to include his own indexable classes in expression math statements, by
// making only reasonable modifications to his own class, without touching
// the xm suite itself.  Here goes...
//===========================================================================//

template<class T>
class UserVec : public xm::Indexable<T,UserVec<T> > {

  public:

    typedef T value_type;
    
  private:
    
    T *v;
    int sz;

    UserVec( const UserVec& );

  public:
    UserVec( int n )
	: xm::Indexable<T,UserVec<T> >(),
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

    T* begin() { return v; }
    T* end() { return v + sz; }

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
 	return assign_from( x );
    }
};

#endif                          // __xm_UserVec_hh__

//---------------------------------------------------------------------------//
//                              end of UserVec.hh
//---------------------------------------------------------------------------//
