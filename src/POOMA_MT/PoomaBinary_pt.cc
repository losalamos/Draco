//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/PoomaBinary_pt.cc
 * \author Randy M. Roberts
 * \date   Thu Nov 18 14:32:28 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Pooma_pt.hh"
#include "PoomaAssign_pt.hh"

#define PETE_TBTREE(OP, T, U) PETE_TBTree<OP, T, U >

#define BinaryAssigner(TYPE, OP) \
Assigner(TYPE, PETE_TBTREE(OP, BFI(TYPE), BFI(TYPE))); \
Assigner(TYPE, PETE_TBTREE(OP, BFI(TYPE), PETE_Scalar<TYPE>)); \
Assigner(VEKTOR(TYPE,6), PETE_TBTREE(OP, BFI(VEKTOR(TYPE,6)), \
                                       PETE_Scalar<TYPE>)); \
Assigner(VEKTOR(TYPE,6), PETE_TBTREE(OP, BFI(VEKTOR(TYPE,6)), \
                                       BFI(VEKTOR(TYPE,6)))); \
Assigner(VEKTOR(TYPE,8), PETE_TBTREE(OP, BFI(VEKTOR(TYPE,8)), \
                                       PETE_Scalar<TYPE>)); \
Assigner(VEKTOR(TYPE,8), PETE_TBTREE(OP, BFI(VEKTOR(TYPE,8)), \
                                       BFI(VEKTOR(TYPE,8))));

BinaryAssigner(double,FnPow);
BinaryAssigner(double,FnArcTan2);
BinaryAssigner(double,FnFmod);

BinaryAssigner(double,OpAdd);
BinaryAssigner(double,OpSubtract);
BinaryAssigner(double,OpDivide);

//---------------------------------------------------------------------------//
//                              end of PoomaBinary_pt.cc
//---------------------------------------------------------------------------//
