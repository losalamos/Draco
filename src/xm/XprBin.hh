//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Created on: Tue Sep 17 09:15:53 1996
// Created by: Geoffrey Furnish
// Also maintained by:
//
//---------------------------------------------------------------------------//

#ifndef __xm_XprBin_hh__
#define __xm_XprBin_hh__

XM_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Binary operations on Indexable's.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class B, class D> \
Xpr< P, XprBinOp<P, ConstRef<P,Indexable<P,A,D> >, ConstRef<P,Indexable<P,B,D> >, ap <P> >, D > \
op ( const Indexable<P,A,D>& a, const Indexable<P,B,D>& b ) \
{ \
    typedef XprBinOp< P, ConstRef<P,Indexable<P,A,D> >, ConstRef<P,Indexable<P,B,D> >, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( ConstRef<P,Indexable<P,A,D> >(a), \
				      ConstRef<P,Indexable<P,B,D> >(b) ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// Junk.
//---------------------------------------------------------------------------//

/*
Keeping this around as an example of how it used to be.  Eventually trash
this. 

template<class P, class A, class B>
Xpr< P, XprBinOp<P, Indexable<P,A>, Indexable<P,B>, OpSub<P> > >
operator-( const Indexable<P,A>& a, const Indexable<P,B>& b )
{
    typedef XprBinOp< P, Indexable<P,A>, Indexable<P,B>, OpSub<P> > ExprT;
    return Xpr< P, ExprT >( ExprT( a, b ) );
}
*/

//---------------------------------------------------------------------------//
// Binary operations between Indexable's and Xpr's.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class B, class D> \
Xpr< P, XprBinOp<P, ConstRef<P,Indexable<P,A,D> >, Xpr<P,B,D>, ap <P> >, D > \
op ( const Indexable<P,A,D>& a, const Xpr<P,B,D>& b ) \
{ \
    typedef XprBinOp< P, ConstRef<P,Indexable<P,A,D> >, Xpr<P,B,D>, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( ConstRef<P,Indexable<P,A,D> >(a), b ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// Binary operations between Xpr's and Indexable's.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class B, class D> \
Xpr< P, XprBinOp<P, Xpr<P,A,D>, ConstRef<P,Indexable<P,B,D> >, ap <P> >, D > \
op ( const Xpr<P,A,D>& a, const Indexable<P,B,D>& b ) \
{ \
    typedef XprBinOp< P, Xpr<P,A,D>, ConstRef<P,Indexable<P,B,D> >, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( a, ConstRef<P,Indexable<P,B,D> >(b) ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// Binary operations between Xpr's and Xpr's.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class B, class D> \
Xpr< P, XprBinOp<P, Xpr<P,A,D>, Xpr<P,B,D>, ap <P> >, D > \
op ( const Xpr<P,A,D>& a, const Xpr<P,B,D>& b ) \
{ \
    typedef XprBinOp< P, Xpr<P,A,D>, Xpr<P,B,D>, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( a, b ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// Val + Idx, and others.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class D> \
Xpr< P, XprBinOp<P, Val<P>, ConstRef<P,Indexable<P,A,D> >, ap <P> >, D > \
op ( const P& a, const Indexable<P,A,D>& b ) \
{ \
    typedef XprBinOp< P, Val<P>, ConstRef<P,Indexable<P,A,D> >, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( Val<P>(a), \
				      ConstRef<P,Indexable<P,A,D> >(b) ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// Idx + Val, and others.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class D> \
Xpr< P, XprBinOp<P, ConstRef<P,Indexable<P,A,D> >, Val<P>, ap <P> >, D > \
op ( const Indexable<P,A,D>& a, const P& b ) \
{ \
    typedef XprBinOp< P, ConstRef<P,Indexable<P,A,D> >, Val<P>, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( ConstRef<P,Indexable<P,A,D> >(a), \
				      Val<P>(b) ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// Val + Xpr, and others.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class D> \
Xpr< P, XprBinOp<P, Val<P>, Xpr<P,A,D>, ap <P> >, D > \
op ( const P& a, const Xpr<P,A,D>& b ) \
{ \
    typedef XprBinOp< P, Val<P>, Xpr<P,A,D>, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( Val<P>(a), b ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// Xpr + Val, and others.
//---------------------------------------------------------------------------//

#define XXX(op,ap) \
template<class P, class A, class D> \
Xpr< P, XprBinOp<P, Xpr<P,A,D>, Val<P>, ap <P> >, D > \
op ( const Xpr<P,A,D>& a, const P& b ) \
{ \
    typedef XprBinOp< P, Xpr<P,A,D>, Val<P>, ap <P> > ExprT; \
    return Xpr< P, ExprT, D >( ExprT( a, Val<P>(b) ) ); \
}

XXX( operator+, OpAdd )
XXX( operator-, OpSub )
XXX( operator*, OpMul )
XXX( operator/, OpDiv )
XXX( pow, OpPow )
XXX( atan2, OpATan2 )
XXX( fmod, OpFmod )
XXX( min, OpMin )
XXX( max, OpMax )

#undef XXX

//---------------------------------------------------------------------------//
// pow.
//---------------------------------------------------------------------------//

template<class P, class A, class D>
Xpr< P, XprBinOp<P, ConstRef<P,Indexable<P,A,D> >, Val<P>, OpPow<P> >, D >
pow( const Indexable<P,A,D>& a, int x )
{
    typedef XprBinOp<P, ConstRef<P,Indexable<P,A,D> >, Val<P>, OpPow<P> > ExprT;
    return Xpr< P, ExprT, D >( ExprT( ConstRef<P,Indexable<P,A,D> >(a),
				      Val<P>( P(x) ) ) );
}

XM_NAMESPACE_END

#endif                          // __xm_XprBin_hh__

//---------------------------------------------------------------------------//
//                              end of xm/XprBin.hh
//---------------------------------------------------------------------------//
