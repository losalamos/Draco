//----------------------------------*-C++-*----------------------------------//
// P13T_PoomaMesh_XYZ_pt.cc
// Randy M. Roberts
// Mon Nov 16 14:18:17 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/P13T.t.cc"
#include "P1Diffusion/SolverP1Diff.hh"
#include "P1Diffusion/P1Diffusion.t.cc"
#include "POOMA_MT/PoomaMesh_XYZ.hh"

typedef UniformCartesian<3> PoomaMesh_t;
typedef PoomaMesh_XYZ<PoomaMesh_t> MT;

typedef rtt_P1Diffusion::SolverP1Diff<MT> MS;

typedef rtt_P1Diffusion::P1Diffusion<MT, MS> DS;

template class rtt_3T::P13T<DS>;

typedef MT::ccsf T1;

#include "timestep/field_ts_advisor.t.cc"

template
void rtt_timestep::field_ts_advisor::update_tstep<T1>(
    const rtt_timestep::ts_manager &, const T1 &, const T1 &);



//---------------------------------------------------------------------------//
//                              end of P13T_PoomaMesh_XYZ_pt.cc
//---------------------------------------------------------------------------//
