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
#include "Global.hh"

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
typename Flat_Mat_State_Builder<MT>::SP_Mat_State
Flat_Mat_State_Builder<MT>::build_Mat_State(SP_Mesh mesh) const
{
    Require (mesh);
    Require (mesh->num_cells() == density.size());

    // make return Mat_State object
    SP_Mat_State return_state;

    // number of cells defined by the mesh 
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar-fields for Mat_State
    typename MT::template CCSF<double> rho(mesh, density);
    typename MT::template CCSF<double> temp(mesh, temperature);
    typename MT::template CCSF<double> sp_heat(mesh, specific_heat);
    
    // create Mat_State object
    return_state = new Mat_State<MT>(rho, temp, sp_heat);

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
 * \param  mesh rtt_dsxx::SP to a mesh
 * \param  mat_state SP to a Mat_State
 * \return SP to an Opacity object
 */
template<class MT>
typename Flat_Mat_State_Builder<MT>::SP_Opacity
Flat_Mat_State_Builder<MT>::build_Opacity(SP_Mesh      mesh,
					  SP_Mat_State mat_state) const
{
    using rtt_mc::global::c;
    using rtt_mc::global::a;

    Require (mesh);
    Require (mat_state);
    Require (mesh->num_cells() == mat_state->num_cells());
    Require (mesh->num_cells() == absorption_opacity.size());

    // return Opacity object
    SP_Opacity return_opacity;
    
    // number of cells
    int num_cells = mesh->num_cells();

    // make cell-centered, scalar fields for Opacity
    typename MT::template CCSF<double> absorption(mesh, absorption_opacity);
    typename MT::template CCSF<double> scattering(mesh, scattering_opacity);
    typename MT::template CCSF<double> fleck(mesh);

    // calculate the Fleck factor in each cell
    double dedT        = 0.0;
    double beta        = 0.0;
    double T           = 0.0;
    double volume      = 0.0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// calculate coefficients needed for Fleck factor
	dedT        = mat_state->get_dedt(cell);
	T           = mat_state->get_T(cell);
	volume      = mesh->volume(cell);
	beta        = 4.0 * a * T*T*T * volume / dedT;
	
	// calculate Fleck factor -- always uses Planck (gray) absorption
	// opacity 
	fleck(cell) = 1.0 / 
	    (1.0 + implicitness * beta * c * delta_t * absorption(cell));
	
	Check (fleck(cell) >= 0.0 && fleck(cell) <= 1.0);
    }
    
    // create Opacity object
    return_opacity = new Opacity<MT>(absorption, scattering, absorption, 
				     fleck);

    Ensure (return_opacity);
    Ensure (return_opacity->num_cells() == num_cells);
    
    return return_opacity;
}

} // end namespace rtt_imc

#endif                          // __imc_Flat_Mat_State_Builder_t_hh__

//---------------------------------------------------------------------------//
//                        end of imc/Flat_Mat_State_Builder.t.hh
//---------------------------------------------------------------------------//
