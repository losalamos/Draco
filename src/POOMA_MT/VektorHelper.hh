//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/VektorHelper.hh
 * \author Randy M. Roberts
 * \date   Tue Nov 16 10:03:40 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __POOMA_MT_VektorHelper_hh__
#define __POOMA_MT_VektorHelper_hh__

#include "AppTypes/Vektor.h"

namespace rtt_POOMA_MT
{

#define VektorHelperUnaryOp(Name, OP)                                        \
template <class T1, unsigned Dim>                                            \
inline Vektor<typename PETEUnaryReturn<T1,OP>::type,Dim>                     \
  Name(const Vektor<T1,Dim> &lhs)                                            \
{                                                                            \
    return TSV_MetaUnary< Vektor<T1,Dim>, OP >::apply(lhs);                  \
}

#define VektorHelperBinaryOp(Name, OP)                                       \
TSV_ELEMENTWISE_OPERATOR(Vektor,Name,OP)

VektorHelperUnaryOp(sin, FnSin);
VektorHelperUnaryOp(cos, FnCos);
VektorHelperUnaryOp(tan, FnTan);
VektorHelperUnaryOp(asin, FnArcSin);
VektorHelperUnaryOp(acos, FnArcCos);
VektorHelperUnaryOp(atan, FnArcTan);
VektorHelperUnaryOp(sinh, FnHypSin);
VektorHelperUnaryOp(cosh, FnHypCos);
VektorHelperUnaryOp(tanh, FnHypTan);
VektorHelperUnaryOp(exp, FnExp);
VektorHelperUnaryOp(log, FnLog);
VektorHelperUnaryOp(log10, FnLog10);
VektorHelperUnaryOp(sqrt, FnSqrt);
VektorHelperUnaryOp(ceil, FnCeil);
VektorHelperUnaryOp(fabs, FnFabs);
VektorHelperUnaryOp(floor, FnFloor);

VektorHelperBinaryOp(pow, FnPow);
VektorHelperBinaryOp(atan2, FnArcTan2);
VektorHelperBinaryOp(fmod, FnFmod);

} // end namespace rtt_POOMA_MT

using rtt_POOMA_MT::sin;
using rtt_POOMA_MT::cos;
using rtt_POOMA_MT::tan;
using rtt_POOMA_MT::asin;
using rtt_POOMA_MT::acos;
using rtt_POOMA_MT::atan;
using rtt_POOMA_MT::sinh;
using rtt_POOMA_MT::cosh;
using rtt_POOMA_MT::tanh;
using rtt_POOMA_MT::exp;
using rtt_POOMA_MT::log;
using rtt_POOMA_MT::log10;
using rtt_POOMA_MT::sqrt;
using rtt_POOMA_MT::ceil;
using rtt_POOMA_MT::fabs;
using rtt_POOMA_MT::floor;

using rtt_POOMA_MT::pow;
using rtt_POOMA_MT::atan2;
using rtt_POOMA_MT::fmod;

// handle abs as a special case.

template<class T1, unsigned Dim>
inline typename PETEUnaryReturn<Vektor<T1, Dim>, FnAbs>::type
PETE_apply(FnAbs, const Vektor<T1, Dim>& a)
{
  return TSV_MetaUnary< Vektor<T1,Dim>, FnAbs>::apply(a);
}

#endif                          // __POOMA_MT_VektorHelper_hh__

//---------------------------------------------------------------------------//
//                              end of POOMA_MT/VektorHelper.hh
//---------------------------------------------------------------------------//
