//----------------------------------*-C++-*----------------------------------//
// MatrixP1Diff_pt.cc
// Randy M. Roberts
// Thu Oct  1 16:06:33 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/MatrixP1Diff.t.cc"
#include "mesh/Mesh_XYZ.hh"
#include "ds++/Mat.hh"

typedef Mesh_XYZ MT;

using namespace rtt_P1Diffusion;

template
class MatrixP1Diff<MT>;

using dsxx::Mat1;

template
void MatrixP1Diff<MT>::multiply(Mat1<double>&, const Mat1<double>&) const;

template
void MatrixP1Diff<MT>::jacobiIteration(Mat1<double>&,
				       const Mat1<double>&) const;

//---------------------------------------------------------------------------//
//                              end of MatrixP1Diff_pt.cc
//---------------------------------------------------------------------------//
