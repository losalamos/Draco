//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Created on: Tue Sep 17 09:16:57 1996
// Created by: Geoffrey Furnish
// Also maintained by:
//
//---------------------------------------------------------------------------//

#ifndef __xm_XprUnary_hh__
#define __xm_XprUnary_hh__

XM_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Provide a variety of operations taking only one argument.
//---------------------------------------------------------------------------//
// Unary operations on Indexable's.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// sin.
//---------------------------------------------------------------------------//

template<class P, class A, class D>
Xpr< P, XprUnaryOp< P, ConstRef<P, Indexable<P,A,D> >, OpSin<P> >, D >
sin( const Indexable<P,A,D>& a )
{
    typedef XprUnaryOp< P, ConstRef<P, Indexable<P,A,D> >, OpSin<P> > ExprT;
    return Xpr< P, ExprT, D >( ExprT( a ) );
}

template<class P, class A, class D>
Xpr< P, XprUnaryOp< P, Xpr<P,A,D>, OpSin<P> >, D >
sin( const Xpr<P,A,D>& a )
{
    typedef XprUnaryOp< P, Xpr<P,A,D>, OpSin<P> > ExprT;
    return Xpr< P, ExprT, D >( ExprT( a ) );
}

XM_NAMESPACE_END

#endif                          // __xm_XprUnary_hh__

//---------------------------------------------------------------------------//
//                              end of xm/XprUnary.hh
//---------------------------------------------------------------------------//
