//----------------------------------*-C++-*----------------------------------//
// MatVecP1Diff_pt.cc
// Randy M. Roberts
// Thu Oct  1 16:03:40 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/MatVecP1Diff.t.cc"
#include "P1Diffusion/MatrixP1Diff.hh"

#include "mesh/Mesh_XYZ.hh"

typedef Mesh_XYZ MT;

using namespace rtt_P1Diffusion;

typedef MatrixP1Diff<MT> Matrix;

template
class MatVecP1Diff<Matrix>;

//---------------------------------------------------------------------------//
//                              end of MatVecP1Diff_pt.cc
//---------------------------------------------------------------------------//
