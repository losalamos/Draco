//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/IMC_Objects_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 22 16:44:10 2003
 * \brief  IMC_Objects_Builder member definitions.
 * \note   Copyright � 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_IMC_Objects_Builder_t_hh
#define rtt_imc_IMC_Objects_Builder_t_hh

#include "ds++/Assert.hh"
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
#include "Hybrid_Diffusion.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
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
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*! 
 * \brief Build the IMC transport objects.
 *
 * This function builds local (on-processor) instantiations of the following
 * objects:
 * - rtt_imc::Gray_Frequency or rtt_imc::Multigroup_Frequency (FT)
 * - rtt_imc::Mat_State
 * - rtt_imc::Opacity
 * - rtt_imc::Diffusion_Opacity (if needed)
 * - rtt_imc::Source
 * - rtt_imc::Tally (with appropriate sub tallies)
 * - rtt_imc::Extrinsic_Surface_Tracker (if needed)
 * - rtt_imc::Random_Walk (if needed)
 * .
 * 
 * The objects should not be inexistence when this function is called.  If it
 * is desired to keep this class in memory use the reset function before
 * calling build_IMC_objects() subsequent times.  Note that reset() does not
 * unset the interface function handed to builder in the constructor.
 *
 * \param mesh rtt_dsxx::SP to a local (on-processor) mesh
 * \param topology rtt_dsxx::SP to a rtt_mc::Topology object
 * \param comm_patterns rtt_dsxx::SP to a rtt_mc::Comm_Patterns object
 * \param rnd_control rtt_dsxx::SP to a rtt_rng::Rnd_Control object
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_IMC_objects(
    SP_Mesh          mesh,
    SP_Topology      topology,
    SP_Comm_Patterns comm_patterns,
    SP_Rnd_Control   rnd_control)
{
    Require (mesh);
    Require (topology);
    Require (comm_patterns);
    Require (rnd_control);

    Require (!frequency);
    Require (!mat_state);
    Require (!opacity);
    Require (!diff_opacity);
    Require (!source);
    Require (!tally);
    Require (!random_walk);
    Require (!tracker);
    Require (!source_builder);

    Require (mesh->num_cells() == topology->num_cells(rtt_c4::node()));

    // first build the material state objects, the appropriate internal
    // function will be called depending upon the base class of the interface
    // (CDI_Data_interface or Flat_Data_Interface)
    build_mat_state_objects(mesh);
    Check (frequency);
    Check (mat_state);
    Check (opacity);

    // build the source
    build_source_object(mesh, topology, rnd_control, comm_patterns);
    Check (source);
    Check (source_builder);

    // build the tally
    build_tally_object(mesh);
    Check (tally);
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
    // build the flat mat state builder
    SP_Mat_State_Builder builder;
    builder = new Flat_Mat_State_Builder<MT,FT>(interface);
    Ensure (builder);
    return builder;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Source<MT,FT,PT> object.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_source_object(
    SP_Mesh          mesh,
    SP_Topology      topology,
    SP_Rnd_Control   rnd_control,
    SP_Comm_Patterns patterns)
{
    Require (mesh);
    Require (topology);
    Require (rnd_control);
    Require (patterns);

    Require (frequency);
    Require (mat_state);
    Require (opacity);

    Require (!source);
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

    // build the source
    source = source_builder->build_Source(mesh, mat_state, opacity,
					  rnd_control, patterns);

    Ensure (source);
    Ensure (source->num_cells() == mesh->num_cells());
    Ensure (source_builder);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Tally<MT> object.
 */
template<class IT, class MT, class FT, class PT>
void IMC_Objects_Builder<IT,MT,FT,PT>::build_tally_object(SP_Mesh mesh)
{
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
	    true);
    Ensure (interface->number_of_surfaces() ? tally->get_Surface_Sub_Tally() :
	    true);
}

} // end namespace rtt_imc

#endif // rtt_imc_IMC_Objects_Builder_t_hh

//---------------------------------------------------------------------------//
//                   end of imc/IMC_Objects_Builder.t.hh
//---------------------------------------------------------------------------//
