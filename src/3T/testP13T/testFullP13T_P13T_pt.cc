#include "3T/testP13T/testFullP13T.hh"
#include "3T/P13T.cc"

typedef Mesh_XYZ MT;
typedef XTM::InterpedMaterialProps MP;
typedef Diffusion_P1<MT> DS;

template class P13T<MT,MP,DS>;
