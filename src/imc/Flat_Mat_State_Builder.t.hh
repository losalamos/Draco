//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Flat_Mat_State_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 14 16:36:43 2001
 * \brief  Flat_Mat_State_Builder member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Flat_Mat_State_Builder_t_hh__
#define __imc_Flat_Mat_State_Builder_t_hh__

#include "Flat_Mat_State_Builder.hh"
#include "Fleck_Factors.hh"
#include "cdi/CDI.hh"
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
	frequency = build_frequency<Dummy_Type>(Type_Switch<FT>());
    }
    Ensure (frequency);

    // build the mat_state
    {
	build_Mat_State(mesh);
    }
    Ensure (mat_state);
    Ensure (mat_state->num_cells() == mesh->num_cells());

    // build the opacity
    {
	opacity = build_opacity<Dummy_Type>(Type_Switch<FT>(), mesh);
    }
    Ensure (opacity);
    Ensure (opacity->num_cells() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build a rtt_imc::Mat_State object.
 *
 * The Mat_State that is returned by this function is defined on the mesh
 * that is input to the function.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 * \return SP to a Mat_State object
 */
template<class MT, class FT>
void Flat_Mat_State_Builder<MT,FT>::build_Mat_State(SP_Mesh mesh)
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
rtt_dsxx::SP<Gray_Frequency> 
Flat_Mat_State_Builder<MT,FT>::build_frequency(Switch_Gray)
{
    // return frequency
    rtt_dsxx::SP<Gray_Frequency> gray_frequency;

    gray_frequency = new Gray_Frequency;

    Ensure (!flat_data->gray_absorption_opacity.empty());
    Ensure (!flat_data->gray_scattering_opacity.empty());
    Ensure (gray_frequency);
    return gray_frequency;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a Multigroup_Frequency.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
rtt_dsxx::SP<Multigroup_Frequency> 
Flat_Mat_State_Builder<MT,FT>::build_frequency(Switch_MG)
{
    Require (!flat_data->group_boundaries.empty());

    // return frequency
    rtt_dsxx::SP<Multigroup_Frequency> mg_frequency;

    mg_frequency = new Multigroup_Frequency(flat_data->group_boundaries);

    Ensure (!flat_data->mg_absorption_opacity.empty());
    Ensure (!flat_data->mg_scattering_opacity.empty());
    Ensure (mg_frequency);
    return mg_frequency;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build an Opacity<MT,Gray_Frequency>.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
rtt_dsxx::SP<Opacity<MT,Gray_Frequency> >
Flat_Mat_State_Builder<MT,FT>::build_opacity(Switch_Gray, SP_Mesh mesh)
{
    using rtt_mc::global::c;
    using rtt_mc::global::a;
    using rtt_dsxx::SP;

    Require (frequency);
    Require (mat_state);
    Require (mesh);

    Check (mesh->num_cells() == flat_data->gray_absorption_opacity.size());
    Check (mesh->num_cells() == flat_data->gray_scattering_opacity.size());

    // return Opacity object
    rtt_dsxx::SP<Opacity<MT,Gray_Frequency> > return_opacity;
    
    // number of cells
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar fields for Opacity
    typename MT::template CCSF<double> absorption(
	mesh, flat_data->gray_absorption_opacity);

    typename MT::template CCSF<double> scattering(
	mesh, flat_data->gray_scattering_opacity);

    // make a Fleck Factors object
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    Check (fleck->fleck.size() == num_cells);

    // calculate the Fleck factor in each cell
    double dedT   = 0.0; // dedT in Jerks/keV
    double T      = 0.0; // temp in keV
    double volume = 0.0; // volume in cc
    double beta   = 0.0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// calculate coefficients needed for Fleck factor
	dedT        = mat_state->get_dedt(cell);
	T           = mat_state->get_T(cell);
	volume      = mesh->volume(cell);

	Check (T      >= 0.0);
	Check (dedT   >  0.0);
	Check (volume >  0.0);

	// calculate beta (4acT^3/Cv)
	beta        = 4.0 * a * T*T*T * volume / dedT;
	
	// calculate Fleck factor
	fleck->fleck(cell) = 1.0 / 
	    (1.0 + implicitness * beta * c * delta_t * absorption(cell));
	
	Check (fleck->fleck(cell) >= 0.0 && fleck->fleck(cell) <= 1.0);
	Check (absorption(cell)   >= 0.0);
	Check (scattering(cell)   >= 0.0);
    }
    
    // create Opacity object
    return_opacity = new Opacity<MT,Gray_Frequency>(frequency, absorption, 
						    scattering, fleck);

    Ensure (return_opacity);
    Ensure (return_opacity->num_cells() == num_cells);
    
    return return_opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build an Opacity<MT,Multigroup_Frequency>.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
rtt_dsxx::SP<Opacity<MT,Multigroup_Frequency> >
Flat_Mat_State_Builder<MT,FT>::build_opacity(Switch_MG, SP_Mesh mesh)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_cdi::CDI;
    using std::pair;
    using rtt_dsxx::soft_equiv;
    using rtt_dsxx::SP;

    Require (frequency);
    Require (mat_state);
    Require (mesh);

    Check (mesh->num_cells() == flat_data->mg_absorption_opacity.size());
    Check (mesh->num_cells() == flat_data->mg_scattering_opacity.size());

    // return opacity
    rtt_dsxx::SP<Opacity<MT,Multigroup_Frequency> > return_opacity;
    
    // number of cells and groups
    int num_cells  = mesh->num_cells();
    int num_groups = frequency->get_num_groups();
    Check (num_groups > 0);

    // make cell-centered, scalar fields for opacities
    typename MT::template CCSF<sf_double> absorption(
	mesh, flat_data->mg_absorption_opacity);
    
    typename MT::template CCSF<sf_double> scattering(
	mesh, flat_data->mg_scattering_opacity);

    typename MT::template CCSF<double>    integrated_norm_planck(mesh);

    typename MT::template CCSF<sf_double> emission_group_cdf(mesh);

    // make a Fleck Factors object
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    Check (fleck->fleck.size() == num_cells);

    // variables needed to calculate the fleck factor and integrated Planck
    // functions 
    double planck = 0.0;
    double dedT   = 0.0;
    double volume = 0.0;
    double beta   = 0.0;
    double T      = 0.0;
    double b_g    = 0.0;

    // loop through cells and build the integrated normalized Planck
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// initialize summation of sigma * b_g over the cell
	double sum = 0.0;

	// resize to the number of groups
	emission_group_cdf(cell).resize(num_groups);

	// get mat state and mesh data
	dedT   = mat_state->get_dedt(cell);
	volume = mesh->volume(cell);
	T      = mat_state->get_T(cell);

	Check (dedT   >  0.0);
	Check (volume >  0.0);
	Check (T      >= 0.0);

	Check (absorption(cell).size() == num_groups);
	Check (scattering(cell).size() == num_groups);

	// calculate the emission group CDF (int sigma * b(x) dx)
	for (int g = 1; g <= num_groups; g++)
	{
	    Check (absorption(cell)[g-1] >= 0.0);
	    Check (scattering(cell)[g-1] >= 0.0);

	    // get the group boundaries
	    pair<double,double> bounds = frequency->get_group_boundaries(g);

	    // integrate the normalized Planckian over the group
	    b_g = CDI::integratePlanckSpectrum(bounds.first, bounds.second,
					       mat_state->get_T(cell)); 

	    // multiply by the absorption opacity and sum
	    sum += b_g * absorption(cell)[g-1];

	    // assign to the cdf
	    emission_group_cdf(cell)[g-1] = sum;
	}

	// integrate the unnormalized Planckian
	integrated_norm_planck(cell) = CDI::integratePlanckSpectrum(
	    frequency->get_group_boundaries().front(),
	    frequency->get_group_boundaries().back(), mat_state->get_T(cell));
	Check (integrated_norm_planck(cell) >= 0.0);

	// calculate the Planckian opacity
	if (integrated_norm_planck(cell) > 0.0)
	    planck = emission_group_cdf(cell).back() /
		integrated_norm_planck(cell); 
	else
	{
	    // weak check that the zero integrated Planck is due to a cold
	    // temperature whose Planckian peak is below the lowest (first)
	    // group boundary.
	    Check (soft_equiv(emission_group_cdf(cell).back(), 0.0));
	    Check (3.0 * mat_state->get_T(cell) 
		   <= frequency->get_group_boundaries(1).first); 

	    // set the ill-defined integrated Planck opacity to zero
	    planck = 0.0;
	}
	Check (planck >= 0.0);

	// calculate beta (4aT^3/Cv)
	beta = 4.0 * a * T*T*T * volume / dedT;
	
	// calculate Fleck Factor
	fleck->fleck(cell) = 1.0 / 
	    (1.0 + implicitness * beta * c * delta_t * planck); 
	Check (fleck->fleck(cell) >= 0.0 && fleck->fleck(cell) <= 1.0);
    }

    // build the return opacity
    return_opacity = new Opacity<MT, Multigroup_Frequency>(
	frequency, absorption, scattering, fleck, integrated_norm_planck, 
	emission_group_cdf);

    Ensure (return_opacity);
    Ensure (return_opacity->num_cells() == mesh->num_cells());
    return return_opacity;
}

} // end namespace rtt_imc

#endif                          // __imc_Flat_Mat_State_Builder_t_hh__

//---------------------------------------------------------------------------//
//                        end of imc/Flat_Mat_State_Builder.t.hh
//---------------------------------------------------------------------------//
