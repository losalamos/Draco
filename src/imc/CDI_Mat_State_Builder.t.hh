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

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR MAT_STATE_BUILDER
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
template<class MT>
typename CDI_Mat_State_Builder<MT>::SP_Mat_State
CDI_Mat_State_Builder<MT>::build_Mat_State(SP_Mesh mesh) const
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    // make return Mat_State object
    SP_Mat_State return_state;

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
    return_state = new Mat_State<MT>(rho, T, Cv);

    Ensure (return_state);
    Ensure (return_state->num_cells() == num_cells);

    // return Mat_State SP
    return return_state;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a rtt_imc::Opacity object.
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
typename CDI_Mat_State_Builder<MT>::SP_Opacity
CDI_Mat_State_Builder<MT>::build_Opacity(SP_Mesh      mesh,
					 SP_Mat_State mat_state) const
{
    using rtt_mc::global::c;
    using rtt_mc::global::a;
    using std::vector;

    Require (mesh);
    Require (mat_state);
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == density.size());

    // return Opacity object
    SP_Opacity return_opacity;
    
    // number of cells
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar fields for Opacity
    typename MT::template CCSF<double> absorption(mesh);
    typename MT::template CCSF<double> scattering(mesh);
    typename MT::template CCSF<double> fleck(mesh);

    // placeholder for CDI objects
    SP_CDI cdi;

    // Reaction types
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;

    // model types
    vector<rtt_cdi::Model> model_abs(material_cdi.size());
    vector<rtt_cdi::Model> model_sct(material_cdi.size());

    // loop through material cdis and determine the model types
    for (int i = 0; i < material_cdi.size(); i++)
    {
	// set SP to cdi
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

    // constants needed for calculation of opacities in /cm
    double T      = 0.0; // temp in keV
    double rho    = 0.0; // density in g/cc
    double dedT   = 0.0; // dedT in Jerks/keV
    double volume = 0.0; // volume in cc
    double beta   = 0.0;

    // cdi index
    int    icdi   = 0;

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

	Check (T      >  0.0);
	Check (rho    >  0.0);
	Check (dedT   >  0.0);
	Check (volume >  0.0);

	// calculate beta (4acT^3/Cv)
	beta   = 4.0 * a * T*T*T * volume / dedT;

	// get the absorption opacity (multiply by density to convert to /cm)
	absorption(cell) = cdi->gray(model_abs[icdi], abs)->
	    getOpacity(T, rho) * rho;

	// get the scattering opacity (multiply by density to convert to /cm)
	scattering(cell) = cdi->gray(model_sct[icdi], sct)->
	    getOpacity(T, rho) * rho;

	// calculate Fleck Factor
	fleck(cell) = 1.0 / 
	    (1.0 + implicitness * beta * c * delta_t * absorption(cell)); 
	
	Check (absorption(cell) >= 0.0);
	Check (scattering(cell) >= 0.0);
	Check (fleck(cell)      >= 0.0 && fleck(cell) <= 1.0);
    }

    // build the return opacity
    return_opacity = new Opacity<MT>(absorption, scattering, fleck);
    
    Ensure (return_opacity);
    Ensure (return_opacity->num_cells() == num_cells);

    return return_opacity;
}

} // end namespace rtt_imc

#endif                          // __imc_CDI_Mat_State_Builder_t_hh__

//---------------------------------------------------------------------------//
//                        end of imc/CDI_Mat_State_Builder.t.hh
//---------------------------------------------------------------------------//
