//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Flat_Mat_State_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 14 16:36:43 2001
 * \brief  Flat_Mat_State_Builder member definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Flat_Mat_State_Builder_t_hh
#define rtt_imc_Flat_Mat_State_Builder_t_hh

#include "Flat_Mat_State_Builder.hh"
#include "Opacity_Builder_Helper.hh"
#include "Fleck_Factors.hh"
#include "ds++/Soft_Equivalence.hh"
#include <utility>

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR MAT_STATE_BUILDER
//---------------------------------------------------------------------------//
/*!
 * \brief Build material state classes.
 *
 * The classes that are built are:
 * - Frequency Type (FT) (either rtt_imc::Gray_Frequency or
 *   rtt_imc::Multigroup_Frequency) 
 * - rtt_imc::Mat_State
 * - rtt_imc::Opacity (gray and multigroup specializations)
 * - rtt_imc::Diffusion_Opacity
 * .
 * \param mesh rtt_dsxx::SP to the mesh
 */
template<class MT, class FT>
void Flat_Mat_State_Builder<MT,FT>::build_mat_classes(SP_Mesh mesh)
{
    using rtt_imc::global::Type_Switch;

    Require (mesh);
    Require (mesh->num_cells() == density.size());
    Require (flat_data);

    Require (!frequency);
    Require (!mat_state);
    Require (!opacity);
    Require (!diff_opacity);

    // build the frequency, specialize on the frequency type
    {
	build_frequency<Dummy_Type>(Type_Switch<FT>());
    }
    Ensure (frequency);

    // build the mat_state
    {
	build_mat_state(mesh);
    }
    Ensure (mat_state);
    Ensure (mat_state->num_cells() == mesh->num_cells());

    // build the opacity
    {
	build_opacity<Dummy_Type>(Type_Switch<FT>(), mesh);
    }
    Ensure (opacity);
    Ensure (opacity->num_cells() == mesh->num_cells());
    Ensure (build_diffusion_opacity ? diff_opacity : true);
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build a rtt_imc::Mat_State object.
 *
 * The Mat_State that is built by this function is defined on the mesh
 * that is input to the function.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 */
template<class MT, class FT>
void Flat_Mat_State_Builder<MT,FT>::build_mat_state(SP_Mesh mesh)
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());
    Check   (mesh->num_cells() == flat_data->specific_heat.size());

    // make cell-centered, scalar-fields for Mat_State
    typename MT::template CCSF<double> rho(mesh, density);
    typename MT::template CCSF<double> temp(mesh, temperature);
    typename MT::template CCSF<double> sp_heat(mesh, flat_data->specific_heat);
    
    // create Mat_State object
    mat_state = new Mat_State<MT>(rho, temp, sp_heat);

    Ensure (mat_state);
    Ensure (mat_state->num_cells() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// PARTIAL SPECIALIZATIONS ON FREQUENCY TYPE (PRIVATE)
//---------------------------------------------------------------------------//
/*!
 * \brief Build a Gray_Frequency.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void Flat_Mat_State_Builder<MT,FT>::build_frequency(Switch_Gray)
{
    // build the gray frequency
    frequency = new Gray_Frequency;

    Ensure (!flat_data->gray_absorption_opacity.empty());
    Ensure (!flat_data->gray_scattering_opacity.empty());
    Ensure (frequency);
    Ensure (frequency->is_gray());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a Multigroup_Frequency.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void Flat_Mat_State_Builder<MT,FT>::build_frequency(Switch_MG)
{
    Require (!flat_data->group_boundaries.empty());

    // build the multigroup frequency
    frequency = new Multigroup_Frequency(flat_data->group_boundaries);

    Ensure (!flat_data->mg_absorption_opacity.empty());
    Ensure (!flat_data->mg_scattering_opacity.empty());
    Ensure (frequency);
    Ensure (frequency->is_multigroup());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build an Opacity<MT,Gray_Frequency>.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void Flat_Mat_State_Builder<MT,FT>::build_opacity(Switch_Gray, SP_Mesh mesh)
{
    using rtt_mc::global::c;
    using rtt_mc::global::a;
    using rtt_dsxx::SP;

    Require (frequency);
    Require (mat_state);
    Require (mesh);

    Check (mesh->num_cells() == flat_data->gray_absorption_opacity.size());
    Check (mesh->num_cells() == flat_data->gray_scattering_opacity.size());
    
    // number of cells
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar fields for Opacity
    typename MT::template CCSF<double> absorption(
	mesh, flat_data->gray_absorption_opacity);

    typename MT::template CCSF<double> scattering(
	mesh, flat_data->gray_scattering_opacity);

    // make a Fleck Factors object
    SP<Fleck_Factors<MT> > fleck =
	Opacity_Builder_Helper<MT,Gray_Frequency>::build_Fleck_Factors(
	    mesh, mat_state, absorption, delta_t, implicitness);
    Check (fleck);
    Check (fleck->fleck.size() == num_cells);
    
    // create Opacity object
    opacity = new Opacity<MT,Gray_Frequency>(frequency, absorption, 
					     scattering, fleck);

    // build the Diffusion_Opacity
    if (build_diffusion_opacity)
    {
	Require (flat_data->rosseland_opacity.size() == num_cells);

	// get the Rosseland opacities
	typename MT::template CCSF<double> rosseland(
	    mesh, flat_data->rosseland_opacity);

	// make the diffusion opacity
	diff_opacity = new Diffusion_Opacity<MT>(fleck, rosseland, scattering);
	Ensure (diff_opacity);
	Ensure (diff_opacity->num_cells() == mesh->num_cells());
    }

    Ensure (opacity);
    Ensure (opacity->num_cells() == num_cells);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build an Opacity<MT,Multigroup_Frequency>.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void Flat_Mat_State_Builder<MT,FT>::build_opacity(Switch_MG, SP_Mesh mesh)
{
    Require (frequency);
    Require (mat_state);
    Require (mesh);

    Check (mesh->num_cells() == flat_data->mg_absorption_opacity.size());
    Check (mesh->num_cells() == flat_data->mg_scattering_opacity.size());

    // make cell-centered, scalar fields for opacities
    typename MT::template CCSF<sf_double> absorption(
	mesh, flat_data->mg_absorption_opacity);
    
    typename MT::template CCSF<sf_double> scattering(
	mesh, flat_data->mg_scattering_opacity);

    // build the opacities
    std::pair<SP_MG_Opacity, SP_Diff_Opacity> opacities = 
	Opacity_Builder_Helper<MT,Multigroup_Frequency>::build_Opacity(
	    mesh, frequency, mat_state, absorption, scattering, delta_t,
	    implicitness, build_diffusion_opacity);
    Check (opacities.first);

    // assign opacities
    opacity      = opacities.first;
    diff_opacity = opacities.second;
    
    Ensure (opacity);
    Ensure (build_diffusion_opacity ? diff_opacity : true);
}

} // end namespace rtt_imc

#endif                          // rtt_imc_Flat_Mat_State_Builder_t_hh

//---------------------------------------------------------------------------//
//                        end of imc/Flat_Mat_State_Builder.t.hh
//---------------------------------------------------------------------------//
