#include "3T/testP13T/testFullP13T.hh"

#include "P1Diffusion/P1Diffusion.t.cc"
#include "P1Diffusion/SolverP1Diff.t.cc"

typedef Mesh_XYZ MT;

using namespace rtt_P1Diffusion;

template
class SolverP1Diff<MT>;

typedef SolverP1Diff<MT> MS;

template
class P1Diffusion<MT,MS>;

#include "P1Diffusion/MatrixP1Diff.t.cc"

template
class MatrixP1Diff<MT>;

typedef MatrixP1Diff<MT> Matrix;

#include "P1Diffusion/MatVecP1Diff.t.cc"

template
class MatVecP1Diff<Matrix>;

#include "P1Diffusion/PreCondP1Diff.t.cc"

template
class PreCondP1Diff<Matrix>;

