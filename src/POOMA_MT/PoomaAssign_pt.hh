//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/PoomaAssign_pt.hh
 * \author Randy M. Roberts
 * \date   Thu Nov 18 14:27:22 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __POOMA_MT_PoomaAssign_pt_hh__
#define __POOMA_MT_PoomaAssign_pt_hh__

#include "Pooma_pt.hh"

#define Assigner(T1, T3) \
template void assign<T1, Dimension_s, T3, OpAssign> \
          (const BareField<T1, Dimension_s> &, T3, OpAssign, \
           ExprTag<true>); \
template void assign<T1, Dimension_s, T3, OpAddAssign> \
          (const BareField<T1, Dimension_s> &, T3, OpAddAssign, \
           ExprTag<true>); \
template void assign<T1, Dimension_s, T3, OpSubtractAssign> \
          (const BareField<T1, Dimension_s> &, T3, OpSubtractAssign, \
           ExprTag<true>); \
template void assign<T1, Dimension_s, T3, OpDivideAssign> \
          (const BareField<T1, Dimension_s> &, T3, OpDivideAssign, \
           ExprTag<true>); \
template void assign<T1, Dimension_s, T3, OpMultiplyAssign> \
          (const BareField<T1, Dimension_s> &, T3, OpMultiplyAssign, \
           ExprTag<true>);

#define BFI(TYPE) BareFieldIterator<TYPE, Dimension_s>

#endif                          // __POOMA_MT_PoomaAssign_pt_hh__

//---------------------------------------------------------------------------//
//                              end of POOMA_MT/PoomaAssign_pt.hh
//---------------------------------------------------------------------------//
