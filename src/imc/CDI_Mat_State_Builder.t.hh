//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/CDI_Mat_State_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Fri Nov 16 11:23:17 2001
 * \brief  CDI_Mat_State_Builder class definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_CDI_Mat_State_Builder_t_hh
#define rtt_imc_CDI_Mat_State_Builder_t_hh

#include "CDI_Mat_State_Builder.hh"
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
void CDI_Mat_State_Builder<MT,FT>::build_mat_classes(SP_Mesh mesh)
{
    using rtt_imc::global::Type_Switch;

    Require (mesh);
    Require (mesh->num_cells() == density.size());

    Require (!frequency);
    Require (!mat_state);
    Require (!opacity);
    Require (!diff_opacity);

    // build the frequency
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

    // build the opacity and diffusion opacity
    {
	build_opacity<Dummy_Type>(Type_Switch<FT>(), mesh);
    }
    Ensure (opacity);
    Ensure (opacity->num_cells() == mesh->num_cells());
    Ensure (build_diffusion_opacity ? diff_opacity : true);
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION FUNCTIONS
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
void CDI_Mat_State_Builder<MT,FT>::build_mat_state(SP_Mesh mesh)
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    // number of cells defined by this mesh
    int num_cells = mesh->num_cells();

    // make CCSF fields required by Mat_State
    typename MT::template CCSF<double> rho(mesh, density);
    typename MT::template CCSF<double> T(mesh, temperature);
    typename MT::template CCSF<double> Cv(mesh);

    // loop through cells and assign specific heats in Jerks/g/keV
    for (int cell = 1; cell <= num_cells; cell++)
    {
	Check (material_cdi[cdi_cell_map[cell-1]-1]->isEoSSet());

	// assign the specific heat from CDI to a temp value and convert from
	// kJ to Jerks
	Cv(cell) = material_cdi[cdi_cell_map[cell-1]-1]->eos()->
	    getElectronHeatCapacity(T(cell), rho(cell)) * 1.0e-6;

	Check (Cv(cell) >= 0.0);
    }

    // create Mat_State object
    mat_state = new rtt_imc::Mat_State<MT>(rho, T, Cv);

    Ensure (mat_state);
    Ensure (mat_state->num_cells() == num_cells);
}

//---------------------------------------------------------------------------//
// PARTIAL SPECIALIZATIONS ON FREQUENCY TYPE (PRIVATE)
//---------------------------------------------------------------------------//
/*!
 * \brief Build a Gray_Frequency.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void CDI_Mat_State_Builder<MT,FT>::build_frequency(Switch_Gray)
{
    // build the gray frequency
    frequency = new Gray_Frequency;

    Ensure (frequency);
    Ensure (frequency->is_gray());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a Multigroup_Frequency.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void CDI_Mat_State_Builder<MT,FT>::build_frequency(Switch_MG)
{
    // build the frequency
    {
	// get the group structure from cdi, ALL multigroup data must have
	// the same group structure, this is enforced by CDI
	sf_double group_bnds = rtt_cdi::CDI::getFrequencyGroupBoundaries();
	frequency            = new Multigroup_Frequency(group_bnds);
    }

    Ensure (frequency);
    Ensure (frequency->get_num_groups() > 0);
    Ensure (frequency->is_multigroup());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a gray rtt_imc::Opacity<MT, Gray_Frequency> object.
 *
 * The Opacity that is built by this function is defined on the mesh that is
 * input to the function.
 *
 * Each CDI opacity contains a rtt_cdi::Model and rtt_cdi::Reaction
 * descriptor i.e., a (rtt_cdi::Model, rtt_cdi::Reaction) pair that is
 * provided by the CDI_Data_Interface.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void CDI_Mat_State_Builder<MT,FT>::build_opacity(Switch_Gray, SP_Mesh mesh)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_dsxx::SP;

    Require (mesh);
    Require (mat_state);
    Require (frequency->is_gray());
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == density.size());
    
    // number of cells
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar fields for Opacity
    typename MT::template CCSF<double> absorption(mesh);
    typename MT::template CCSF<double> scattering(mesh);

    // constants needed for calculation of opacities in /cm
    double T      = 0.0; // temp in keV
    double rho    = 0.0; // density in g/cc

    // cdi index
    int    icdi   = 0;
    SP_CDI cdi;

    // reaction types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;

    // models
    rtt_cdi::Model model_abs;
    rtt_cdi::Model model_sct;

    // loop through the cells and assign the opacities and fleck factor
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// cdi index
	icdi = cdi_cell_map[cell-1]-1;

	// set SP to cdi
	cdi = material_cdi[icdi];
	Check (cdi);

	// get material data
	T      = mat_state->get_T(cell);
	rho    = mat_state->get_rho(cell);
	Check (T   >=  0.0);
	Check (rho >   0.0);

	// get the models (the standard allows explicit casts from int to a
	// valid enumeration)
	model_abs = static_cast<rtt_cdi::Model>(cdi_models[icdi].first);
	model_sct = static_cast<rtt_cdi::Model>(cdi_models[icdi].second);

	Check(model_abs == rtt_cdi::ROSSELAND ||
	      model_abs == rtt_cdi::PLANCK    ||
	      model_abs == rtt_cdi::ANALYTIC);
	Check(model_sct == rtt_cdi::ROSSELAND ||
	      model_sct == rtt_cdi::PLANCK    ||
	      model_sct == rtt_cdi::ANALYTIC);

	Check (cdi->isGrayOpacitySet(model_abs, abs));
	Check (cdi->isGrayOpacitySet(model_sct, sct));

	// get the absorption opacity (multiply by density to convert to /cm)
	absorption(cell) = cdi->gray(model_abs, abs)->getOpacity(T, rho) * rho;

	// get the scattering opacity (multiply by density to convert to /cm)
	scattering(cell) = cdi->gray(model_sct, sct)->getOpacity(T, rho) * rho; 
	
	Check (absorption(cell)   >= 0.0);
	Check (scattering(cell)   >= 0.0);
    }

    // build fleck factors
    SP_Fleck_Factors fleck =
	Opacity_Builder_Helper<MT,Gray_Frequency>::build_Fleck_Factors(
	    mesh, mat_state, absorption, delta_t, implicitness);
    Check (fleck);
    Check (fleck->fleck.size() == num_cells);

    // build the return opacity
    opacity = new Opacity<MT, Gray_Frequency>(frequency, absorption, 
					      scattering, fleck);

    // build the diffusion opacity if required
    if (build_diffusion_opacity)
    {
	// build the Rosseland opacities
	build_diff_opacity_gray<Dummy_Type>(Switch_Gray(), mesh, fleck);
	Ensure (diff_opacity);
	Ensure (diff_opacity->num_cells() == mesh->num_cells());
    }

    Ensure (opacity);
    Ensure (opacity->num_cells() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Build a Diffusion_opacity object for gray problems.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void CDI_Mat_State_Builder<MT,FT>::build_diff_opacity_gray(
    Switch_Gray,
    SP_Mesh          mesh,
    SP_Fleck_Factors fleck)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_dsxx::SP;

    Require (frequency->is_gray());
    Require (mesh);
    Require (mat_state);
    Require (fleck);
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == density.size());
    
    // number of cells
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar fields for Rosseland opacity
    typename MT::template CCSF<double> rosseland(mesh);

    // constants needed for calculation of opacities in /cm
    double T      = 0.0; // temp in keV
    double rho    = 0.0; // density in g/cc

    // cdi index
    int    icdi   = 0;
    SP_CDI cdi;

    // reaction and model types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Model    ros = rtt_cdi::ROSSELAND;

    // loop through the cells and assign the opacities and fleck factor
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// cdi index
	icdi = cdi_cell_map[cell-1]-1;

	// set SP to cdi
	cdi = material_cdi[icdi];
	Check (cdi);

	// get material data
	T   = mat_state->get_T(cell);
	rho = mat_state->get_rho(cell);

	Check (T   >= 0.0);
	Check (rho >  0.0);

	// the Rosseland absorption opacity must be defined in the CDI
	Check (cdi->isGrayOpacitySet(ros, abs));
	
	// get the rosseland absorption from cdi
	// (multiply by density to convert to /cm)
	rosseland(cell) = cdi->gray(ros, abs)->getOpacity(T, rho) * rho;

	Check (rosseland(cell) >= 0.0);
    }

    // build the diffusion opacity
    diff_opacity = new Diffusion_Opacity<MT>(fleck, rosseland);
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Build a mulitgroup rtt_imc::Opacity<MT, Multigroup_Frequency>
 * object.
 *
 * The Opacity that is built by this function is defined on the mesh that is
 * input to the function.
 *
 * Each CDI opacity contains a rtt_cdi::Model and rtt_cdi::Reaction
 * descriptor i.e., a (rtt_cdi::Model, rtt_cdi::Reaction) pair.  The CDI
 * objects are queried to determine the model type for the absorption and
 * scattering opacities.  Even though CDI allows multiple models for each
 * reaction type, the CDI_Mat_State_Builder expects CDI to use only one model
 * for each reaction type.  For example, CDI_Mat_State_Builder will throw an
 * exception if it encounters opacities defined by (rtt_cdi::PLANCK,
 * rtt_cdi::ABSORPTION) and (rtt_cdi::ANALYTIC, rtt_cdi::ABSORPTION) in the
 * same CDI object.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void CDI_Mat_State_Builder<MT,FT>::build_opacity(Switch_MG, SP_Mesh mesh)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_cdi::CDI;
    using std::pair;
    using rtt_dsxx::soft_equiv;
    using rtt_dsxx::SP;

    Require (mesh);
    Require (mat_state);
    Require (frequency);
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == density.size());
    Require (CDI::getNumberFrequencyGroups() == frequency->get_num_groups());
 
    // make cell-centered, scalar fields for opacities
    typename MT::template CCSF<sf_double> absorption(mesh);
    typename MT::template CCSF<sf_double> scattering(mesh);

    // get the multigroup absorption and scattering opacities in /cm,
    // integrate the Planckian
    build_mg_opacities<Dummy_Type>(Switch_MG(), absorption, scattering);

    // build the opacities
    std::pair<SP_Opacity, SP_Diff_Opacity> opacities = 
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

//---------------------------------------------------------------------------//
/*!
 * \brief Get multigroup absorption and scattering opacities.
 */
template<class MT, class FT>
template<class Stop_Explicit_Instantiation>
void CDI_Mat_State_Builder<MT,FT>::build_mg_opacities(Switch_MG,
						      ccvf &absorption,
						      ccvf &scattering)
{
    Require (mat_state);
    Require (frequency->is_multigroup());

    // get the number of cells
    int num_cells  = absorption.get_Mesh().num_cells();
    int num_groups = rtt_cdi::CDI::getNumberFrequencyGroups();

    Check (mat_state->num_cells() == num_cells);

    // constants needed for calculation of opacities in /cm
    double T   = 0.0; // temp in keV
    double rho = 0.0; // density in g/cc

    // cdi index
    int    icdi = 0;
    SP_CDI cdi;

    // reaction types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;

    // models
    rtt_cdi::Model model_abs;
    rtt_cdi::Model model_sct;

    // loop through the cells and assign the opacities 
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// cdi index
	icdi = cdi_cell_map[cell-1]-1;

	// set SP to cdi
	cdi = material_cdi[icdi];
	Check (cdi);

	// get material data
	T   = mat_state->get_T(cell);
	rho = mat_state->get_rho(cell);

	Check (T   >= 0.0);
	Check (rho >  0.0);

	// get the models (the standard allows explicit casts from int to a
	// valid enumeration)
	model_abs = static_cast<rtt_cdi::Model>(cdi_models[icdi].first);
	model_sct = static_cast<rtt_cdi::Model>(cdi_models[icdi].second);

	Check(model_abs == rtt_cdi::ROSSELAND ||
	      model_abs == rtt_cdi::PLANCK    ||
	      model_abs == rtt_cdi::ANALYTIC);
	Check(model_sct == rtt_cdi::ROSSELAND ||
	      model_sct == rtt_cdi::PLANCK    ||
	      model_sct == rtt_cdi::ANALYTIC);
	
	// get the absorption coefficients in cm^2/g 
	absorption(cell) = cdi->mg(model_abs, abs)->getOpacity(T, rho);

	// get the scattering coefficients in cm^2/g 
	scattering(cell) = cdi->mg(model_sct, sct)->getOpacity(T, rho);

	Check (absorption(cell).size() == num_groups);
	Check (scattering(cell).size() == num_groups);
	
	// loop through groups and convert to /cm by multiplying by rho
	for (int i = 0; i < num_groups; i++)
	{
	    absorption(cell)[i] *= rho;
	    scattering(cell)[i] *= rho;

	    Check (absorption(cell)[i] >= 0.0);
	    Check (scattering(cell)[i] >= 0.0);
	}
    }
}

} // end namespace rtt_imc

#endif                          // rtt_imc_CDI_Mat_State_Builder_t_hh

//---------------------------------------------------------------------------//
//                        end of imc/CDI_Mat_State_Builder.t.hh
//---------------------------------------------------------------------------//
