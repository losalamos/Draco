//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

#ifndef __xm_Xpr_hh__
#define __xm_Xpr_hh__

#ifndef __xm_xm_hh__
#error "Users should only include xm/xm.hh"
#endif

XM_NAMESPACE_BEG

//===========================================================================//
// class ConstRef - 

// 
//===========================================================================//

template<class P, class T>
class ConstRef {
    const T& t;
  public:
    ConstRef( const T& _t ) : t(_t) {}
    P operator[]( int n ) const { return t[n]; }
};

//===========================================================================//
// class Val - 

// 
//===========================================================================//

template<class P>
class Val {
    P v;
  public:
    explicit Val( const P& _v ) : v(_v) {}

    P operator[]( int n ) const { return v; }
};

//===========================================================================//
// class Xpr - 

// 
//===========================================================================//

template<class P, class E, class D>
class Xpr
{
    E e;
  public:
    Xpr( const E& _e ) : e(_e) {}

    P operator[]( int n ) const { return e[n]; }
};

//===========================================================================//
// class Indexable - 

// 
//===========================================================================//

template<class P, class I, class D = I>
class Indexable {
  public:
    explicit Indexable() {}
    Indexable( const Indexable& a ) {}

    P operator[]( int n ) const
    {
	return static_cast<const I*>(this)->operator[](n);
    }

    template<class E>
    I& assign_from( const Xpr<P,E,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) = x[i];

 	return *me;
    }

// Operators taking a Val.
    Indexable& operator+=( const P& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) += x;

	return *this;
    }

    Indexable& operator-=( const P& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) -= x;

	return *this;
    }

    Indexable& operator*=( const P& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) *= x;

	return *this;
    }

    Indexable& operator/=( const P& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) /= x;

	return *this;
    }

// Operators taking an Indexable.
    Indexable& operator+=( const Indexable<P,I,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) += x[i];

	return *this;
    }

    Indexable& operator-=( const Indexable<P,I,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) -= x[i];

	return *this;
    }

    Indexable& operator*=( const Indexable<P,I,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) *= x[i];

	return *this;
    }

    Indexable& operator/=( const Indexable<P,I,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) /= x[i];

	return *this;
    }

// Operators taking an Xpr.
    template<class E>
    Indexable& operator+=( const Xpr<P,E,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) += x[i];

	return *this;
    }

    template<class E>
    Indexable& operator-=( const Xpr<P,E,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) -= x[i];

	return *this;
    }

    template<class E>
    Indexable& operator*=( const Xpr<P,E,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) *= x[i];

	return *this;
    }

    template<class E>
    Indexable& operator/=( const Xpr<P,E,D>& x )
    {
	I *me = static_cast<I*>(this);

	for( int i=0; i < me->size(); i++ )
	    me->operator[](i) /= x[i];

	return *this;
    }
};

//===========================================================================//
// class XprBinOp - 

// 
//===========================================================================//

template<class P, class A, class B, class Op>
class XprBinOp {
    A a;
    B b;
  public:
    XprBinOp( const A& _a, const B& _b )
	: a(_a), b(_b) {}

    P operator[]( int n ) const
    {
	return Op::apply( a[n], b[n] );
    }
};

//===========================================================================//
// class XprUnaryOp - 

// 
//===========================================================================//

template<class P, class A, class Op>
class XprUnaryOp {
    A a;
  public:
    XprUnaryOp( const A& _a ) : a(_a) {}

    P operator[]( int n ) const
    {
	return Op::apply( a[n] );
    }
};

XM_NAMESPACE_END

#include "xm/XprBin.hh"
#include "xm/XprUnary.hh"

#endif                          // __xm_Xpr_hh__

//---------------------------------------------------------------------------//
//                              end of xm/Xpr.hh
//---------------------------------------------------------------------------//
