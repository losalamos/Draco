//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/CDI_Mat_State_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Fri Nov 16 11:23:17 2001
 * \brief  CDI_Mat_State_Builder class definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_CDI_Mat_State_Builder_t_hh__
#define __imc_CDI_Mat_State_Builder_t_hh__

#include "CDI_Mat_State_Builder.hh"
#include "Frequency.hh"
#include "Opacity.hh"
#include "Mat_State.hh"
#include "Global.hh"
#include "Fleck_Factors.hh"

namespace rtt_imc
{

//===========================================================================//
/*!
 * \class CDI_Mat_State_Builder_Helper
 *
 * \brief Helps CDI_Mat_State_Builder build Mat_State objects.
 *
 * This is a helper class used by specializations of CDI_Mat_State_Builder to
 * build Mat_State objects.  We use this class because specializations are
 * not required to build the Mat_State; thus, each class specialization of
 * Mat_State_Builder can use the same implementation.
 *
 * The CDI_Mat_State_Builder class
 */
//===========================================================================//

template<class MT>
class CDI_Mat_State_Builder_Helper
{
    public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>             SP_Mesh;
    typedef rtt_dsxx::SP<rtt_cdi::CDI>   SP_CDI;
    typedef std::vector<SP_CDI>          sf_CDI;
    typedef std::vector<int>             sf_int;
    typedef std::vector<double>          sf_double;    
    typedef rtt_dsxx::SP<Mat_State<MT> > SP_Mat_State;
    
    // Unbound friendship declaration.
    template<class M, class F>
    friend class CDI_Mat_State_Builder;

  private:
    // Build a Mat_State.
    static SP_Mat_State build_Mat_State(SP_Mesh, const sf_CDI &,
					const sf_int &, const sf_double &,
					const sf_double &);
};

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Mat_State class for all specializations of <MT,FT> types.
 */
template<class MT>
typename CDI_Mat_State_Builder_Helper<MT>::SP_Mat_State 
CDI_Mat_State_Builder_Helper<MT>::build_Mat_State(
    SP_Mesh          mesh,
    const sf_CDI    &material_cdi,
    const sf_int    &cdi_cell_map,
    const sf_double &density,
    const sf_double &temperature)
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    // make return Mat_State object
    rtt_dsxx::SP<rtt_imc::Mat_State<MT> > return_state;

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
    return_state = new rtt_imc::Mat_State<MT>(rho, T, Cv);

    Ensure (return_state);
    Ensure (return_state->num_cells() == num_cells);

    // return Mat_State SP
    return return_state;
}

//===========================================================================//
// CDI_MAT_STATE_BUILDER<MT,GRAY_FREQUENCY>
//===========================================================================//
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR MAT_STATE_BUILDER
//---------------------------------------------------------------------------//
/*!
 * \brief Build a rtt_imc::Gray_Frequency object.
 */
template<class MT>
typename CDI_Mat_State_Builder<MT,Gray_Frequency>::SP_Frequency
CDI_Mat_State_Builder<MT,Gray_Frequency>::build_Frequency()
{
    Require (material_cdi.size() > 0);

    // return Frequency type
    SP_Frequency frequency(new Gray_Frequency);

    Ensure (frequency);
    Ensure (material_cdi[0]->isGrayOpacitySet(
		rtt_cdi::PLANCK, rtt_cdi::ABSORPTION) ||
	    material_cdi[0]->isGrayOpacitySet(
		rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION) ||
	    material_cdi[0]->isGrayOpacitySet(
		rtt_cdi::ANALYTIC, rtt_cdi::ABSORPTION));

    return frequency;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a rtt_imc::Mat_State object for CDI_Mat_State_Builder<MT,
 * Gray_Frequency> specialization.
 *
 * The Mat_State that is returned by this function is defined on the mesh
 * that is input to the function.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 * \return SP to a Mat_State object
 */
template<class MT>
typename CDI_Mat_State_Builder<MT,Gray_Frequency>::SP_Mat_State
CDI_Mat_State_Builder<MT, Gray_Frequency>::build_Mat_State(
    SP_Mesh mesh)
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    // make return Mat_State object
    SP_Mat_State return_state = CDI_Mat_State_Builder_Helper<MT>::
	build_Mat_State(mesh, material_cdi, cdi_cell_map, density,
			temperature);

    Ensure (return_state);
    Ensure (return_state->num_cells() == mesh->num_cells());
    return return_state;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a gray rtt_imc::Opacity<MT, Gray_Frequency> object.
 *
 * The Opacity that is returned by this function is defined on the mesh that
 * is input to the function.
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
 * \param  mat_state SP to a Mat_State
 * \return SP to an Opacity object
 */
template<class MT>
typename CDI_Mat_State_Builder<MT,Gray_Frequency>::SP_Opacity
CDI_Mat_State_Builder<MT,Gray_Frequency>::build_Opacity(
    SP_Mesh      mesh,
    SP_Frequency freq,
    SP_Mat_State mat_state)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_dsxx::SP;

    Require (mesh);
    Require (mat_state);
    Require (freq);
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == density.size());

    // return opacity
    SP_Opacity return_opacity;
    
    // number of cells
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar fields for Opacity
    typename MT::template CCSF<double> absorption(mesh);
    typename MT::template CCSF<double> scattering(mesh);

    // make a fleck factors object
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    Check (fleck->fleck.size() == num_cells);

    // determine the model types for gray opacity
    sf_Opacity_Model model_abs(material_cdi.size());
    sf_Opacity_Model model_sct(material_cdi.size());
    determine_gray_models(model_abs, model_sct);

    // constants needed for calculation of opacities in /cm
    double T      = 0.0; // temp in keV
    double rho    = 0.0; // density in g/cc
    double dedT   = 0.0; // dedT in Jerks/keV
    double volume = 0.0; // volume in cc
    double beta   = 0.0;

    // cdi index
    int    icdi   = 0;
    SP_CDI cdi;

    // reaction types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;

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
	dedT   = mat_state->get_dedt(cell);
	volume = mesh->volume(cell);

	Check (T      >=  0.0);
	Check (rho    >   0.0);
	Check (dedT   >   0.0);
	Check (volume >   0.0);

	// calculate beta (4acT^3/Cv)
	beta   = 4.0 * a * T*T*T * volume / dedT;

	// get the absorption opacity (multiply by density to convert to /cm)
	absorption(cell) = cdi->gray(model_abs[icdi], abs)->
	    getOpacity(T, rho) * rho;

	// get the scattering opacity (multiply by density to convert to /cm)
	scattering(cell) = cdi->gray(model_sct[icdi], sct)->
	    getOpacity(T, rho) * rho;

	// calculate Fleck Factor
	fleck->fleck(cell) = 1.0 / 
	    (1.0 + implicitness * beta * c * delta_t * absorption(cell)); 
	
	Check (absorption(cell)   >= 0.0);
	Check (scattering(cell)   >= 0.0);
	Check (fleck->fleck(cell) >= 0.0 && fleck->fleck(cell) <= 1.0);
    }

    // build the return opacity
    return_opacity = new Opacity<MT, Gray_Frequency>(
	freq, absorption, scattering, fleck);

    Ensure (return_opacity);
    Ensure (return_opacity->num_cells() == mesh->num_cells());
    return return_opacity;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Detemine the models used for gray opacities.
 */
template<class MT>
void CDI_Mat_State_Builder<MT, Gray_Frequency>::determine_gray_models(
    sf_Opacity_Model &model_abs, 
    sf_Opacity_Model &model_sct)
{
    Require (model_abs.size() == material_cdi.size());
    Require (model_sct.size() == material_cdi.size());
    Require (material_cdi.size() > 0);

    // placeholder for CDI objects
    SP_CDI cdi;

    // Reaction types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;

    // loop over each material and determine the (model, reaction) combos for
    // each 
    for (int i = 0; i < material_cdi.size(); i++)
    {
	// get a SP to the cdi for this material
	cdi = material_cdi[i];
	Check (cdi);

	// absorption opacities, we check to make sure only ONE of the
	// following (model, ABSORPTION) combinations is defined

	if (cdi->isGrayOpacitySet(rtt_cdi::ANALYTIC, abs))
	{
	    // use Analytic values for the absorption and Planck opacities
	    model_abs[i]    = rtt_cdi::ANALYTIC;

	    // make sure the other absorption models are not set
	    Insist (!cdi->isGrayOpacitySet(rtt_cdi::ROSSELAND, abs) &&
		    !cdi->isGrayOpacitySet(rtt_cdi::PLANCK, abs),
		    "Tried to define multiple opacity models for absorption.");
	}
	else if (cdi->isGrayOpacitySet(rtt_cdi::ROSSELAND, abs))
	{
	    // the absorption opacity is Rosseland
	    model_abs[i]    = rtt_cdi::ROSSELAND;

	    // make sure the other absorption models are not set
	    Insist (!cdi->isGrayOpacitySet(rtt_cdi::PLANCK, abs),
		    "Tried to define multiple opacity models for absorption.");
	}
	else if (cdi->isGrayOpacitySet(rtt_cdi::PLANCK, abs))
	{
	    // the absorption opacity is Planck;
	    model_abs[i]    = rtt_cdi::PLANCK;
	}
	else
	{
	    Insist (0, "No absorption opacity available!");
	}

	// scattering opacities, we check to make sure only ONE of the
	// following (model, SCATTERING) combinations is defined

	if (cdi->isGrayOpacitySet(rtt_cdi::ANALYTIC, sct))
	{
	    model_sct[i] = rtt_cdi::ANALYTIC;

	    // make sure the other scattering models are not set
	    Insist (!cdi->isGrayOpacitySet(rtt_cdi::ROSSELAND, sct) &&
		    !cdi->isGrayOpacitySet(rtt_cdi::PLANCK, sct),
		    "Tried to define multiple opacity models for scattering.");
	}
	else if (cdi->isGrayOpacitySet(rtt_cdi::ROSSELAND, sct))
	{
	    model_sct[i] = rtt_cdi::ROSSELAND;

	    // make sure the other absorption models are not set
	    Insist (!cdi->isGrayOpacitySet(rtt_cdi::PLANCK, abs),
		    "Tried to define multiple opacity models for scattering.");
	}
	else if (cdi->isGrayOpacitySet(rtt_cdi::PLANCK, sct))
	{
	    model_sct[i] = rtt_cdi::PLANCK;
	}
	else
	{
	    Insist (0, "No scattering opacity available.");
	}
    }
}

//===========================================================================//
// CDI_MAT_STATE_BUILDER<MT,MULTIGROUP_FREQUENCY>
//===========================================================================//
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR MAT_STATE_BUILDER
//---------------------------------------------------------------------------//
/*!
 * \brief Build a rtt_imc::Multigroup_Frequency object.
 */
template<class MT>
typename CDI_Mat_State_Builder<MT,Multigroup_Frequency>::SP_Frequency
CDI_Mat_State_Builder<MT,Multigroup_Frequency>::build_Frequency()
{
    Require (material_cdi.size() > 0);
    Require (material_cdi[0]->isMultigroupOpacitySet(
		 rtt_cdi::PLANCK, rtt_cdi::ABSORPTION) ||
	     material_cdi[0]->isMultigroupOpacitySet(
		 rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION) ||
	     material_cdi[0]->isMultigroupOpacitySet(
		 rtt_cdi::ANALYTIC, rtt_cdi::ABSORPTION));

    // get the group structure from cdi, ALL multigroup data must have the
    // same group structure, this is enforced by CDI
    sf_double group_bnds = rtt_cdi::CDI::getFrequencyGroupBoundaries();

    // return Frequency type
    SP_Frequency frequency(new Multigroup_Frequency(group_bnds));

    Ensure (frequency);
    Ensure (frequency->get_num_groups() > 0);
    return frequency;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a rtt_imc::Mat_State object for CDI_Mat_State_Builder<MT,
 * Multigroup_Frequency> specialization.
 *
 * The Mat_State that is returned by this function is defined on the mesh
 * that is input to the function.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 * \return SP to a Mat_State object
 */
template<class MT>
typename CDI_Mat_State_Builder<MT,Multigroup_Frequency>::SP_Mat_State
CDI_Mat_State_Builder<MT, Multigroup_Frequency>::build_Mat_State(
    SP_Mesh mesh)
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    // make return Mat_State object
    SP_Mat_State return_state = CDI_Mat_State_Builder_Helper<MT>::
	build_Mat_State(mesh, material_cdi, cdi_cell_map, density,
			temperature);

    Ensure (return_state);
    Ensure (return_state->num_cells() == mesh->num_cells());
    return return_state;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Build a mulitgroup rtt_imc::Opacity<MT, Multigroup_Frequency>
 * object.
 *
 * The Opacity that is returned by this function is defined on the mesh that
 * is input to the function.
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
 * \param  mat_state SP to a Mat_State
 * \return SP to an Opacity object
 */
template<class MT>
typename CDI_Mat_State_Builder<MT,Multigroup_Frequency>::SP_Opacity
CDI_Mat_State_Builder<MT,Multigroup_Frequency>::build_Opacity(
    SP_Mesh      mesh,
    SP_Frequency freq,
    SP_Mat_State mat_state)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_cdi::CDI;
    using rtt_dsxx::soft_equiv;
    using rtt_dsxx::SP;

    Require (mesh);
    Require (mat_state);
    Require (freq);
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == density.size());
    Require (CDI::getNumberFrequencyGroups() == freq->get_num_groups());

    // return opacity
    SP_Opacity return_opacity;
    
    // number of cells and groups
    int num_cells  = mesh->num_cells();
    int num_groups = freq->get_num_groups();
 
    // make cell-centered, scalar fields for opacities
    typename MT::template CCSF<sf_double> absorption(mesh);
    typename MT::template CCSF<sf_double> scattering(mesh);
    typename MT::template CCSF<double>    integrated_norm_planck(mesh);
    typename MT::template CCSF<sf_double> emission_group_cdf(mesh);

    // make a Fleck Factors object
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    Check (fleck->fleck.size() == num_cells);

    // get the multigroup absorption and scattering opacities in /cm,
    // integrate the Planckian
    build_opacities(mat_state, absorption, scattering);

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

	// calculate the emission group CDF (int sigma * b(x) dx)
	for (int g = 1; g <= num_groups; g++)
	{
	    // integrate the normalized Planckian over the group
	    b_g = CDI::integratePlanckSpectrum(g, mat_state->get_T(cell)); 

	    // multiply by the absorption opacity and sum
	    sum += b_g * absorption(cell)[g-1];

	    // assign to the cdf
	    emission_group_cdf(cell)[g-1] = sum;
	}

	// integrate the unnormalized Planckian
	integrated_norm_planck(cell) = 
	    CDI::integratePlanckSpectrum(mat_state->get_T(cell));
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
		   <= freq->get_group_boundaries(1).first); 

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
	freq, absorption, scattering, fleck, integrated_norm_planck, 
	emission_group_cdf);

    Ensure (return_opacity);
    Ensure (return_opacity->num_cells() == mesh->num_cells());
    return return_opacity;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Determine the models used for multigroup opacities.
 */
template<class MT>
void CDI_Mat_State_Builder<MT,Multigroup_Frequency>::determine_mg_models(
    sf_Opacity_Model &model_abs, 
    sf_Opacity_Model &model_sct)
{
    Require (model_abs.size() == material_cdi.size());
    Require (model_sct.size() == material_cdi.size());
    Require (material_cdi.size() > 0);

    // placeholder for CDI objects
    SP_CDI cdi;

    // Reaction types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;

    // loop over each material and determine the (model, reaction) combos for
    // each 
    for (int i = 0; i < material_cdi.size(); i++)
    {
	// get a SP to the cdi for this material
	cdi = material_cdi[i];
	Check (cdi);

	// absorption opacities, we check to make sure only ONE of the
	// following (model, ABSORPTION) combinations is defined

	if (cdi->isMultigroupOpacitySet(rtt_cdi::ANALYTIC, abs))
	{
	    // use Analytic values for the absorption and Planck opacities
	    model_abs[i]    = rtt_cdi::ANALYTIC;

	    // make sure the other absorption models are not set
	    Insist (!cdi->isMultigroupOpacitySet(rtt_cdi::ROSSELAND, abs) &&
		    !cdi->isMultigroupOpacitySet(rtt_cdi::PLANCK, abs),
		    "Tried to define multiple opacity models for absorption.");
	}
	else if (cdi->isMultigroupOpacitySet(rtt_cdi::ROSSELAND, abs))
	{
	    // the absorption opacity is Rosseland
	    model_abs[i]    = rtt_cdi::ROSSELAND;

	    // make sure the other absorption models are not set
	    Insist (!cdi->isMultigroupOpacitySet(rtt_cdi::PLANCK, abs),
		    "Tried to define multiple opacity models for absorption.");
	}
	else if (cdi->isMultigroupOpacitySet(rtt_cdi::PLANCK, abs))
	{
	    // the absorption opacity is Planck;
	    model_abs[i]    = rtt_cdi::PLANCK;
	}
	else
	{
	    Insist (0, "No absorption opacity available!");
	}

	// scattering opacities, we check to make sure only ONE of the
	// following (model, SCATTERING) combinations is defined

	if (cdi->isMultigroupOpacitySet(rtt_cdi::ANALYTIC, sct))
	{
	    model_sct[i] = rtt_cdi::ANALYTIC;

	    // make sure the other scattering models are not set
	    Insist (!cdi->isMultigroupOpacitySet(rtt_cdi::ROSSELAND, sct) &&
		    !cdi->isMultigroupOpacitySet(rtt_cdi::PLANCK, sct),
		    "Tried to define multiple opacity models for scattering.");
	}
	else if (cdi->isMultigroupOpacitySet(rtt_cdi::ROSSELAND, sct))
	{
	    model_sct[i] = rtt_cdi::ROSSELAND;

	    // make sure the other absorption models are not set
	    Insist (!cdi->isMultigroupOpacitySet(rtt_cdi::PLANCK, abs),
		    "Tried to define multiple opacity models for scattering.");
	}
	else if (cdi->isMultigroupOpacitySet(rtt_cdi::PLANCK, sct))
	{
	    model_sct[i] = rtt_cdi::PLANCK;
	}
	else
	{
	    Insist (0, "No scattering opacity available.");
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get multigroup absorption and scattering opacities.
 */
template<class MT>
void CDI_Mat_State_Builder<MT,Multigroup_Frequency>::build_opacities(
    SP_Mat_State                           mat_state,
    typename MT::template CCSF<sf_double> &absorption,
    typename MT::template CCSF<sf_double> &scattering)
{
    // get the number of cells
    int num_cells  = absorption.get_Mesh().num_cells();
    int num_groups = rtt_cdi::CDI::getNumberFrequencyGroups();

    Check (mat_state->num_cells() == num_cells);

    // determine the model types for multigroup opacity
    sf_Opacity_Model model_abs(material_cdi.size());
    sf_Opacity_Model model_sct(material_cdi.size());
    determine_mg_models(model_abs, model_sct);

    // constants needed for calculation of opacities in /cm
    double T   = 0.0; // temp in keV
    double rho = 0.0; // density in g/cc

    // cdi index
    int    icdi   = 0;
    SP_CDI cdi;

    // reaction types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;

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
	
	// get the absorption coefficients in cm^2/g 
	absorption(cell) = cdi->mg(model_abs[icdi], abs)->getOpacity(T, rho);

	// get the scattering coefficients in cm^2/g 
	scattering(cell) = cdi->mg(model_sct[icdi], sct)->getOpacity(T, rho);

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

#endif                          // __imc_CDI_Mat_State_Builder_t_hh__

//---------------------------------------------------------------------------//
//                        end of imc/CDI_Mat_State_Builder.t.hh
//---------------------------------------------------------------------------//
