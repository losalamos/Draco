//----------------------------------*-C++-*----------------------------------//
// Parallel_Builder_pt.cc
// Thomas M. Evans
// Tue Apr 14 15:31:12 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Parallel_Builder
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "Particle.hh"
#include "Parallel_Builder.t.hh"

IMCSPACE

using RNG::Rnd_Control;

typedef Particle_Buffer<Particle<OS_Mesh> > POS_Buffer;
typedef Source<OS_Mesh> OS_Source;

template class Parallel_Builder<OS_Mesh>;

template
void Parallel_Builder<OS_Mesh>::dist_census(const Source_Init<OS_Mesh> &, 
					    const POS_Buffer &,
					    POS_Buffer::Census &);

template
void Parallel_Builder<OS_Mesh>::recv_census(const POS_Buffer &,
					    POS_Buffer::Census &);

template SP<OS_Source> 
Parallel_Builder<OS_Mesh>::send_Source(SP<OS_Mesh>, SP<Mat_State<OS_Mesh> >, 
				       SP<Rnd_Control>, 
				       const Source_Init<OS_Mesh> &, 
				       const POS_Buffer &);

template SP<OS_Source > 
Parallel_Builder<OS_Mesh>::recv_Source(SP<OS_Mesh>, SP<Mat_State<OS_Mesh> >, 
				       SP<Rnd_Control>, const POS_Buffer &);

template SP<Communicator<Particle<OS_Mesh> > > 
Parallel_Builder<OS_Mesh>::build_Communicator(int);

template SP<Communicator<Particle<OS_Mesh> > > 
Parallel_Builder<OS_Mesh>::send_Communicator();

template SP<Communicator<Particle<OS_Mesh> > > 
Parallel_Builder<OS_Mesh>::recv_Communicator();

CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Builder_pt.cc
//---------------------------------------------------------------------------//
