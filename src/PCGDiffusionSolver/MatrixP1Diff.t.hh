//----------------------------------*-C++-*----------------------------------//
// MatrixP1Diff.t.cc
// Randy M. Roberts
// Mon Sep 28 08:33:48 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "MatrixP1Diff.hh"
#include "ds++/Assert.hh"

namespace rtt_PCGDiffusionSolver
{

 template<class MT>
 void MatrixP1Diff<MT>::multiply(ccsf &b, const ccsf &x) const
 {
     b = diagonal_m * x;

     MT::gather(xFC, x, MT::OpAssign());
     MT::swap_faces(xSwap, xFC, 0.0);

     bFC = offDiagonal_m * xSwap;

     MT::scatter(b, bFC, MT::OpAddAssign());
     
 }
 
 template<class MT>
 template<class FT>
 void MatrixP1Diff<MT>::multiply(FT &b, const FT &x) const
 {
     Assert(btmp.size() == b.size());
     Assert(xtmp.size() == x.size());

     std::copy(x.begin(), x.end(), xtmp.begin());

     multiply(btmp, xtmp);

     std::copy(btmp.begin(), btmp.end(), b.begin());
 }
 
 template<class MT>
 void MatrixP1Diff<MT>::jacobiIteration(ccsf &b, const ccsf &x) const
 {
     b = x / diagonal_m;
 }
 
 template<class MT>
 template<class FT>
 void MatrixP1Diff<MT>::jacobiIteration(FT &b, const FT &x) const
 {
     Assert(btmp.size() == b.size());
     Assert(xtmp.size() == x.size());

     std::copy(x.begin(), x.end(), xtmp.begin());

     jacobiIteration(btmp, xtmp);

     std::copy(btmp.begin(), btmp.end(), b.begin());
 }
 
} // end namespace rtt_PCGDiffusionSolver


//---------------------------------------------------------------------------//
//                              end of MatrixP1Diff.t.cc
//---------------------------------------------------------------------------//
