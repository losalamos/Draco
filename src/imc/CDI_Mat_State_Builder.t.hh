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
#include "Global.hh"
#include "ds++/Soft_Equivalence.hh"
#include <utility>

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
 * \brief Build material state classes.
 *
 * This specialization builds material state data for gray problems.  The
 * classes that are built are:
 * - rtt_imc::Gray_Frequency
 * - rtt_imc::Mat_State
 * - rtt_imc::Opacity (gray specialization)
 * - rtt_imc::Diffusion_Opacity
 * .
 * \param mesh rtt_dsxx::SP to the mesh
 */
template<class MT>
void CDI_Mat_State_Builder<MT,Gray_Frequency>::build_mat_classes(SP_Mesh mesh)
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    Require (!gray);
    Require (!mat_state);
    Require (!opacity);
    Require (!diff_opacity);

    // build the frequency
    {
	gray = new Gray_Frequency;
    }
    Ensure (gray);

    // build the mat_state
    {
	mat_state = CDI_Mat_State_Builder_Helper<MT>::build_Mat_State(
	    mesh, material_cdi, cdi_cell_map, density, temperature);
    }
    Ensure (mat_state);
    Ensure (mat_state->num_cells() == mesh->num_cells());

    // build the opacity and diffusion opacity
    {
	build_Opacity(mesh);
    }
    Ensure (opacity);
    Ensure (opacity->num_cells() == mesh->num_cells());
    Ensure (build_diffusion_opacity ? diff_opacity : true);
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build a gray rtt_imc::Opacity<MT, Gray_Frequency> object.
 *
 * The Opacity that is returned by this function is defined on the mesh that
 * is input to the function.
 *
 * Each CDI opacity contains a rtt_cdi::Model and rtt_cdi::Reaction
 * descriptor i.e., a (rtt_cdi::Model, rtt_cdi::Reaction) pair that is
 * provided by the CDI_Data_Interface.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 */
template<class MT>
void CDI_Mat_State_Builder<MT,Gray_Frequency>::build_Opacity(SP_Mesh mesh)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_dsxx::SP;

    Require (mesh);
    Require (mat_state);
    Require (gray);
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
    SP<Fleck_Factors<MT> > fleck = 
	Opacity_Builder_Helper<MT,Gray_Frequency>::build_Fleck_Factors(
	    mesh, mat_state, absorption, delta_t, implicitness);
    Check (fleck);
    Check (fleck->fleck.size() == num_cells);

    // build the return opacity
    opacity = new Opacity<MT, Gray_Frequency>(gray, absorption, scattering, 
					      fleck);

    // build the diffusion opacity if required
    if (build_diffusion_opacity)
    {
	// build the Rosseland opacities
	build_Diffusion_Opacity(mesh, fleck);
	Ensure (diff_opacity);
	Ensure (diff_opacity->num_cells() == mesh->num_cells());
    }

    Ensure (opacity);
    Ensure (opacity->num_cells() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Build a Diffusion_opacity object.
 *
 * \param  mesh rtt_dsxx::SP to a mesh
 */
template<class MT>
void CDI_Mat_State_Builder<MT,Gray_Frequency>::build_Diffusion_Opacity(
    SP_Mesh                          mesh,
    rtt_dsxx::SP<Fleck_Factors<MT> > fleck)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_dsxx::SP;

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
    rtt_cdi::Reaction tot = rtt_cdi::TOTAL;
    rtt_cdi::Reaction sct = rtt_cdi::SCATTERING;
    rtt_cdi::Reaction abs = rtt_cdi::ABSORPTION;
    rtt_cdi::Model    ros = rtt_cdi::ROSSELAND;
    rtt_cdi::Model    mabs;
    rtt_cdi::Model    msct;

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

	// determine how to calculate the Rosseland opacity:
	// (multiply by density to convert to /cm)
	
	// if the Rosseland total opacity is given then use it
	if (cdi->isGrayOpacitySet(ros, tot))
	{
	    rosseland(cell) = cdi->gray(ros, tot)->getOpacity(T, rho) * rho;
	}

	// otherwise add the absorption and scattering together to estimate
	// the Rosseland opacity (assuming that the absorption is not Planck)
	else if (cdi_models[icdi].first != rtt_cdi::PLANCK)
	{
	    // get the models for the absorption and scattering
	    mabs = static_cast<rtt_cdi::Model>(cdi_models[icdi].first);
	    msct = static_cast<rtt_cdi::Model>(cdi_models[icdi].second);
	    
	    // now add scattering to absorption to get the Rosseland total
	    rosseland(cell) = (cdi->gray(mabs, abs)->getOpacity(T, rho)  +
			       cdi->gray(msct, sct)->getOpacity(T, rho)) * rho;
	}

	// if the only thing we are given is the Planck opacity for
	// absorption then we fail
	else
	{
	    throw rtt_dsxx::assertion(
		"Cannot build Rosseland opacity from Planck absorption.");
	}

	Check (rosseland(cell) >= 0.0);
    }

    // build the diffusion opacity
    diff_opacity = new Diffusion_Opacity<MT>(fleck, rosseland);
}

//===========================================================================//
// CDI_MAT_STATE_BUILDER<MT,MULTIGROUP_FREQUENCY>
//===========================================================================//
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR MAT_STATE_BUILDER
//---------------------------------------------------------------------------//
/*!
 * \brief Build material state classes.
 *
 * This specialization builds material state data for multigroup problems.
 * The classes that are built are:
 * - rtt_imc::Multigroup_Frequency
 * - rtt_imc::Mat_State
 * - rtt_imc::Opacity (gray specialization)
 * - rtt_imc::Diffusion_Opacity
 * .
 * \param mesh rtt_dsxx::SP to the mesh
 */
template<class MT>
void CDI_Mat_State_Builder<MT,Multigroup_Frequency>::build_mat_classes(
    SP_Mesh mesh)
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    Require (!mg);
    Require (!mat_state);
    Require (!opacity);
    Require (!diff_opacity);

    // build the frequency
    {
	// get the group structure from cdi, ALL multigroup data must have
	// the same group structure, this is enforced by CDI
	sf_double group_bnds = rtt_cdi::CDI::getFrequencyGroupBoundaries();
	mg                   = new Multigroup_Frequency(group_bnds);
    }
    Ensure (mg);
    Ensure (mg->get_num_groups() > 0);

    // build the mat_state
    {
	mat_state = CDI_Mat_State_Builder_Helper<MT>::build_Mat_State(
	    mesh, material_cdi, cdi_cell_map, density, temperature);
    }
    Ensure (mat_state);
    Ensure (mat_state->num_cells() == mesh->num_cells());

    // build the opacity
    {
	build_Opacity(mesh);
    }
    Ensure (opacity);
    Ensure (opacity->num_cells() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION FUNCTIONS
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
void CDI_Mat_State_Builder<MT,Multigroup_Frequency>::build_Opacity(
    SP_Mesh mesh)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_cdi::CDI;
    using std::pair;
    using rtt_dsxx::soft_equiv;
    using rtt_dsxx::SP;

    Require (mesh);
    Require (mat_state);
    Require (mg);
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == density.size());
    Require (CDI::getNumberFrequencyGroups() == mg->get_num_groups());
 
    // make cell-centered, scalar fields for opacities
    typename MT::template CCSF<sf_double> absorption(mesh);
    typename MT::template CCSF<sf_double> scattering(mesh);

    // get the multigroup absorption and scattering opacities in /cm,
    // integrate the Planckian
    build_opacities(absorption, scattering);

    // build the opacities
    std::pair<SP_Opacity, SP_Diff_Opacity> opacities = 
	Opacity_Builder_Helper<MT,Multigroup_Frequency>::build_Opacity(
	    mesh, mg, mat_state, absorption, scattering, delta_t,
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
template<class MT>
void CDI_Mat_State_Builder<MT,Multigroup_Frequency>::build_opacities(
    typename MT::template CCSF<sf_double> &absorption,
    typename MT::template CCSF<sf_double> &scattering)
{
    Require (mat_state);

    // get the number of cells
    int num_cells  = absorption.get_Mesh().num_cells();
    int num_groups = rtt_cdi::CDI::getNumberFrequencyGroups();

    Check (mat_state->num_cells() == num_cells);

    // constants needed for calculation of opacities in /cm
    double T   = 0.0; // temp in keV
    double rho = 0.0; // density in g/cc

    // cdi index
    int    icdi   = 0;
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
