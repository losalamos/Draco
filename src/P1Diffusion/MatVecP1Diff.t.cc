//----------------------------------*-C++-*----------------------------------//
// MatVecP1Diff.t.ccc
// Randy M. Roberts
// Fri Sep 25 13:43:12 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/MatVecP1Diff.hh"

#include <iostream>

namespace rtt_P1Diffusion {

 template<class Matrix>
 void MatVecP1Diff<Matrix>::MatVec(Mat1 &b, const Mat1 &x)
 {
     spMatrix->multiply(b, x);
 }

} // end namespace

//---------------------------------------------------------------------------//
//                              end of MatVecP1Diff.t.cc
//---------------------------------------------------------------------------//
