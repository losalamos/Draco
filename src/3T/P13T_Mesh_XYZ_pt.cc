//----------------------------------*-C++-*----------------------------------//
// P13T_Mesh_XYZ_pt.cc
// Randy M. Roberts
// Tue Nov  3 14:56:33 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/P13T.t.cc"
#include "P1Diffusion/SolverP1Diff.hh"
#include "P1Diffusion/P1Diffusion.hh"
#include "mesh/Mesh_XYZ.hh"

typedef Mesh_XYZ MT;

typedef rtt_P1Diffusion::SolverP1Diff<MT> MS;

typedef rtt_P1Diffusion::P1Diffusion<MT, MS> DS;

template class rtt_3T::P13T<DS>;

typedef MT::ccsf T1;

#include "timestep/field_ts_advisor.t.cc"

template
void rtt_timestep::field_ts_advisor::update_tstep<T1>(
    const rtt_timestep::ts_manager &, const T1 &, const T1 &);

//---------------------------------------------------------------------------//
//                              end of P13T_Mesh_XYZ_pt.cc
//---------------------------------------------------------------------------//
