//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/IMC_Objects_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 22 16:44:10 2003
 * \brief  IMC_Objects_Builder member definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_IMC_Objects_Builder_t_hh
#define rtt_imc_IMC_Objects_Builder_t_hh

#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "rng/Random.hh"
#include "mc/Comm_Patterns.hh"
#include "mc/Topology.hh"
#include "CDI_Mat_State_Builder.hh"
#include "Flat_Mat_State_Builder.hh"
#include "CDI_Data_Interface.hh"
#include "Flat_Data_Interface.hh"
#include "IMC_Objects_Builder.hh"
#include "Rep_Source_Builder.hh"
#include "DD_Source_Builder.hh"
#include "Tally_Builder.hh"
#include "Extrinsic_Tracker_Builder.hh"
#include "Hybrid_Diffusion.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR AND DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param iface interface that is a derived class of
 * rtt_imc::CDI_Data_Interface or rtt_imc::Flat_Data_Interface,
 * rtt_imc::Interface, and rtt_imc::Surface_Tracking_Interface
 */
template<class IT, class MT, class FT, class PT>
IMC_Objects_Builder<IT,MT,FT,PT>::IMC_Objects_Builder(SP_Interface iface)
    : interface(iface)
{
    Require (interface);
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build the material state objects.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_mat_state_objects(SP_Mesh mesh)
{
    using rtt_dsxx::SP;

    Require (interface);
    Require (mesh);

    Require (!frequency);
    Require (!mat_state);
    Require (!opacity);
    Require (!diff_opacity);

    // get a mat state builder, we get the appropriate builder through the
    // interface base class (CDI_Mat_State_Builder if interface isA
    // CDI_Data_Interface or Flat_Mat_State_Builder if interface isA
    // Flat_Data_Interface)
    SP<Mat_State_Builder<MT,FT> > builder = get_mat_builder<IT>(*interface);
    Check (builder);

    // build the mat objects
    builder->build_mat_classes(mesh);

    // assign objects

    // get the frequency
    frequency    = builder->get_Frequency();
	
    // get the mat_state
    mat_state    = builder->get_Mat_State();

    // get the opacity 
    opacity      = builder->get_Opacity();

    // get the diffusion opacity
    diff_opacity = builder->get_Diffusion_Opacity();

    // if we need hybrid diffusion then diff opacity better exist
    Ensure (interface->get_hybrid_diffusion_method() ? diff_opacity : 
	    !diff_opacity);

    Ensure (frequency);
    Ensure (mat_state);
    Ensure (opacity);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a CDI_Mat_State_Builder.
 *
 * We add an extra template argument here because we cannot instantiate this
 * explicitly if the interface does not provide CDI_Data_Interface.
 */
template<class IT, class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
typename IMC_Objects_Builder<IT,MT,FT,PT>::SP_Mat_State_Builder
IMC_Objects_Builder<IT,MT,FT,PT>::get_mat_builder(
    const CDI_Data_Interface &i)
{
    Require (interface);

    // build the cdi mat state builder
    SP_Mat_State_Builder builder;
    builder = new CDI_Mat_State_Builder<MT,FT>(interface);
    Ensure (builder);
    return builder;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a Flat_Mat_State_Builder.
 *
 * We add an extra template argument here because we cannot instantiate this
 * explicitly if the interface does not provide Flat_Data_Interface.
 */
template<class IT, class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
typename IMC_Objects_Builder<IT,MT,FT,PT>::SP_Mat_State_Builder
IMC_Objects_Builder<IT,MT,FT,PT>::get_mat_builder(
    const Flat_Data_Interface &i)
{
    Require (interface);

    // build the flat mat state builder
    SP_Mat_State_Builder builder;
    builder = new Flat_Mat_State_Builder<MT,FT>(interface);
    Ensure (builder);
    return builder;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Source_Builder<MT,FT,PT> object.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_source_builder_object(
    SP_Mesh          mesh,
    SP_Topology      topology)
{
    Require (interface);
    Require (mesh);
    Require (topology);

    Require (frequency);
    Require (mat_state);
    Require (opacity);

    Require (!source_builder);

    // build the source builder based on the topology
    if (topology->get_parallel_scheme() == "replication")
    {
	// point source_builder to a Rep_Source_Builder
	source_builder =
	    new Rep_Source_Builder<MT,FT,PT>(interface, mesh, topology);
    }
    else if (topology->get_parallel_scheme() == "DD")
    {
	// point source_builder to a DD_Source_Builder
	source_builder = 
	    new DD_Source_Builder<MT,FT,PT>(interface, mesh, topology);
    }
    else
    {
	// haven't added any other sources yet
	Insist(0, "Don't have other source support yet!");
    }

    Ensure (source_builder);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Source<MT,FT,PT> object.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_source_object(
    SP_Mesh          mesh,
    SP_Rnd_Control   rnd_control,
    SP_Comm_Patterns patterns)
{
    Require (interface);
    Require (mesh);
    Require (rnd_control);
    Require (patterns);

    Require (frequency);
    Require (mat_state);
    Require (opacity);

    Require (source_builder);

    Require (!source);

    // build the source
    source = source_builder->build_Source(mesh, mat_state, opacity,
					  rnd_control, patterns);

    Ensure (source);
    Ensure (source->num_cells() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Tally<MT> object.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_tally_object(SP_Mesh mesh)
{
    Require (interface);
    Require (mesh);

    Require (!tally);

    // make a tally builder
    Tally_Builder<MT> builder(interface);
    
    // build the tally
    tally = builder.build_Tally(mesh);

    Ensure (tally);

    // check to make sure that the tally has the appropriate sub tallies
    // defined 
    Ensure (interface->get_hybrid_diffusion_method() 
	    == Hybrid_Diffusion::RANDOM_WALK ? tally->get_RW_Sub_Tally() :
	    !tally->get_RW_Sub_Tally());
    Ensure (interface->number_of_surfaces() ? tally->get_Surface_Sub_Tally() :
	    !tally->get_Surface_Sub_Tally());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build time-dependent objects that determine particle transport
 * behavior.
 *
 * Objects built by this function must be built every cycle.  The objects
 * built by this function are:
 * - rtt_imc::Random_Walk
 * .
 * The interface determines whether these should be built or not.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_time_dependent_particle_objects(
    SP_Mesh mesh)
{
    Require (interface);
    Require (mesh);

    Require (!random_walk);

    // build the random walk
    if (interface->get_hybrid_diffusion_method() 
	== Hybrid_Diffusion::RANDOM_WALK)
    {
	Require (diff_opacity);
	random_walk = new Random_Walk<MT>(mesh, diff_opacity);
    }

    Ensure (interface->get_hybrid_diffusion_method() 
	    == Hybrid_Diffusion::RANDOM_WALK ? random_walk : !random_walk);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build time-independent objects that determine particle transport
 * behavior.
 *
 * Objects built by this function do not necessarily have to be built every
 * cycle.  However, they may be if problem options change. The objects built
 * by this function are:
 * - rtt_imc::Extrinsic_Surface_Tracker
 * .
 * The interface determines whether these should be built or not.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_time_independent_particle_objects(
    SP_Mesh mesh)
{
    Require (interface);
    Require (mesh);

    Require (!tracker);

    // build the extrinsic surface tracker
    {
	Extrinsic_Tracker_Builder<MT> builder(*mesh, interface);
	tracker = builder.build_tracker();
    }

    // the extrinsic surface tracker may not exist on a decomposed part of
    // the mesh even though surface tracking is on; thus we can't easily
    // check this but we'll try
    Remember (int count = 0;);
    Remember (if (tracker) count = 1;);
    Remember (rtt_c4::global_prod(count););
    Ensure   (interface->number_of_surfaces() ? count : !count);
}

} // end namespace rtt_imc

#endif // rtt_imc_IMC_Objects_Builder_t_hh

//---------------------------------------------------------------------------//
//                   end of imc/IMC_Objects_Builder.t.hh
//---------------------------------------------------------------------------//
