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
// Unary operations on Indexable's.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class D> \
Xpr< P, XprUnaryOp< P, ConstRef<P, Indexable<P,A,D> >, ap <P> >, D > \
op ( const Indexable<P,A,D>& a ) \
{ \
    typedef XprUnaryOp< P, ConstRef<P, Indexable<P,A,D> >, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( a ) ); \
}

XXX( operator+, OpUnaryAdd )
XXX( operator-, OpUnarySub )
XXX( sin, OpSin )
XXX( cos, OpCos )
XXX( tan, OpTan )
XXX( asin, OpASin )
XXX( acos, OpACos )
XXX( atan, OpATan )
XXX( cosh, OpCosh )
XXX( sinh, OpSinh )
XXX( tanh, OpTanh )
XXX( exp, OpExp )
XXX( log, OpLog )
XXX( log10, OpLog10 )
XXX( sqrt, OpSqrt )
XXX( ceil, OpCeil )
XXX( abs, OpAbs )
XXX( labs, OpLabs )
XXX( fabs, OpFabs )
XXX( floor, OpFloor )

#undef XXX

//---------------------------------------------------------------------------//
// Unary operations on Xpr's.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class D> \
Xpr< P, XprUnaryOp< P, Xpr<P,A,D>, ap <P> >, D > \
op ( const Xpr<P,A,D>& a ) \
{ \
    typedef XprUnaryOp< P, Xpr<P,A,D>, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( a ) ); \
}

XXX( operator+, OpUnaryAdd )
XXX( operator-, OpUnarySub )
XXX( sin, OpSin )
XXX( cos, OpCos )
XXX( tan, OpTan )
XXX( asin, OpASin )
XXX( acos, OpACos )
XXX( atan, OpATan )
XXX( cosh, OpCosh )
XXX( sinh, OpSinh )
XXX( tanh, OpTanh )
XXX( exp, OpExp )
XXX( log, OpLog )
XXX( log10, OpLog10 )
XXX( sqrt, OpSqrt )
XXX( ceil, OpCeil )
XXX( abs, OpAbs )
XXX( labs, OpLabs )
XXX( fabs, OpFabs )
XXX( floor, OpFloor )

#undef XXX

XM_NAMESPACE_END

#endif                          // __xm_XprUnary_hh__

//---------------------------------------------------------------------------//
//                              end of xm/XprUnary.hh
//---------------------------------------------------------------------------//
