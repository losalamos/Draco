//----------------------------------*-C++-*----------------------------------//
// MatrixP1Diff.t.cc
// Randy M. Roberts
// Mon Sep 28 08:33:48 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "MatrixP1Diff.hh"
#include "ds++/Assert.hh"

namespace rtt_P1Diffusion
{

 template<class MT>
 void MatrixP1Diff<MT>::multiply(ccsf &b, const ccsf &x) const
 {
     b = *spDiagonal * x;

     MT::gather(*spxFC, x, MT::OpAssign());
     MT::swap_faces(*spxSwap, *spxFC);

     *spbFC = *spOffDiagonal * (*spxSwap);

     MT::scatter(b, *spbFC, MT::OpAddAssign());
 }
 
 template<class MT>
 template<class FT>
 void MatrixP1Diff<MT>::multiply(FT &b, const FT &x) const
 {
     Assert(spbtmp->size() == b.size());
     Assert(spxtmp->size() == x.size());

     std::copy(x.begin(), x.end(), spxtmp->begin());

     multiply(*spbtmp, *spxtmp);

     std::copy(spbtmp->begin(), spbtmp->end(), b.begin());
 }
 
 template<class MT>
 void MatrixP1Diff<MT>::jacobiIteration(ccsf &b, const ccsf &x) const
 {
     b = x / *spDiagonal;
 }
 
 template<class MT>
 template<class FT>
 void MatrixP1Diff<MT>::jacobiIteration(FT &b, const FT &x) const
 {
     Assert(spbtmp->size() == b.size());
     Assert(spxtmp->size() == x.size());

     std::copy(x.begin(), x.end(), spxtmp->begin());

     jacobiIteration(*spbtmp, *spxtmp);

     std::copy(spbtmp->begin(), spbtmp->end(), b.begin());
 }
 
} // end namespace rtt_P1Diffusion


//---------------------------------------------------------------------------//
//                              end of MatrixP1Diff.t.cc
//---------------------------------------------------------------------------//
