//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/PoomaUnary_pt.cc
 * \author Randy M. Roberts
 * \date   Thu Nov 18 14:29:01 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Pooma_pt.hh"
#include "PoomaAssign_pt.hh"

#define PETE_TUTREE(OP, TYPE) PETE_TUTree<OP, TYPE >

#define UnaryVektorAssigner(TYPE, SIZE, OP) \
Assigner(VEKTOR(TYPE,SIZE), PETE_TUTREE(OP, BFI(VEKTOR(TYPE,SIZE))))

#define UnaryAssigner(TYPE, OP) \
Assigner(TYPE, PETE_TUTREE(OP, BFI(TYPE))) \
UnaryVektorAssigner(TYPE, 6, OP); \
UnaryVektorAssigner(TYPE, 8, OP);


UnaryAssigner(double,FnSin);
UnaryAssigner(double,FnCos);
UnaryAssigner(double,FnTan);
UnaryAssigner(double,FnArcSin);
UnaryAssigner(double,FnArcCos);
UnaryAssigner(double,FnArcTan);
UnaryAssigner(double,FnHypSin);
UnaryAssigner(double,FnHypCos);
UnaryAssigner(double,FnHypTan);
UnaryAssigner(double,FnExp);
UnaryAssigner(double,FnLog);
UnaryAssigner(double,FnLog10);
UnaryAssigner(double,FnSqrt);
UnaryAssigner(double,FnCeil);
UnaryAssigner(double,FnFabs);
UnaryAssigner(double,FnFloor);

UnaryAssigner(double,OpUnaryPlus);
UnaryAssigner(double,OpUnaryMinus);

UnaryAssigner(int,FnAbs);

//---------------------------------------------------------------------------//
//                              end of PoomaUnary_pt.cc
//---------------------------------------------------------------------------//
