//----------------------------------*-C++-*----------------------------------//
// PreCondP1Diff_pt.cc
// Randy M. Roberts
// Thu Oct  1 16:05:32 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "PreCondP1Diff.t.hh"
#include "MatrixP1Diff.hh"
#include "mesh/Mesh_XYZ.hh"

typedef Mesh_XYZ MT;

using namespace rtt_P1Diffusion;

typedef MatrixP1Diff<MT> Matrix;

template
class PreCondP1Diff<Matrix>;

//---------------------------------------------------------------------------//
//                              end of PreCondP1Diff_pt.cc
//---------------------------------------------------------------------------//
