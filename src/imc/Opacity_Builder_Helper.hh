//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Opacity_Builder_Helper.hh
 * \author Thomas M. Evans
 * \date   Fri Aug  8 16:18:04 2003
 * \brief  Opacity_Builder_Helper class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Opacity_Builder_Helper_hh
#define rtt_imc_Opacity_Builder_Helper_hh

#include "Frequency.hh"
#include "Opacity.hh"
#include "Fleck_Factors.hh"
#include "Diffusion_Opacity.hh"
#include "Mat_State.hh"
#include "Global.hh"
#include "cdi/CDI.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <utility>
#include <vector>

namespace rtt_imc
{

//===========================================================================//
/*!
 * \class Opacity_Builder_Helper
 * \brief Helps Mat_State_Builder classes build opacity data.
 *
 * This class is used by the Mat_State_Builder hierarchy
 * (CDI_Mat_State_Builder and Flat_Mat_State_Builder) to build opacity data.
 * In effect, it takes functionality that is common to both derived builder
 * classes and implements it.  We do this using a helper class instead of the
 * Mat_State_Builder base class because it is easier to specialize on FT
 * (Frequency Type) with a helper class.
 *
 * For gray specializations the helper builds Fleck factors.  For multigroup
 * specializations the helper builds rtt_imc::Opacity and
 * rtt_imc::Diffusion_Opacity objects.
 */
/*!
 * \example imc/test/tstOpacity_Builder_Helper.cc
 *
 * Test of Opacity_Builder_Helper class.
 */
// revision history:
// -----------------
// 0) (Fri Aug  8 16:18:04 2003) Thomas M. Evans: original
// 1) 25-AUG-2003 : updated to integrate the Rosseland absorption opacity
//                  instead of the Rosseland total opacity
// 
//===========================================================================//

template<class MT, class FT>
class Opacity_Builder_Helper 
{
  private:
    Opacity_Builder_Helper();
};

//===========================================================================//
/*!
 * \class Opacity_Builder_Helper<MT,Gray_Frequency>
 *
 * \brief Specialization of Opacity_Builder_Helper on Gray_Frequency.
 */
//===========================================================================//

template<class MT>
class Opacity_Builder_Helper<MT, Gray_Frequency>
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>                   SP_Mesh;
    typedef rtt_dsxx::SP<Fleck_Factors<MT> >   SP_Fleck;
    typedef rtt_dsxx::SP<Mat_State<MT> >       SP_Mat_State;
    typedef typename MT::template CCSF<double> ccsf;

  public:
    // Build Fleck Factors.
    static SP_Fleck build_Fleck_Factors(SP_Mesh, SP_Mat_State, const ccsf &,
					const double, const double);
};

//---------------------------------------------------------------------------//
/*!
 * \brief Build Fleck Factors for Opacity<MT,Gray> creation.
 */
template<class MT>
typename Opacity_Builder_Helper<MT,Gray_Frequency>::SP_Fleck 
Opacity_Builder_Helper<MT,Gray_Frequency>::build_Fleck_Factors(
    SP_Mesh       mesh,
    SP_Mat_State  mat_state,
    const ccsf   &absorption,
    const double  delta_t,
    const double  implicitness)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    
    Require (mesh);
    Require (mat_state);

    // number of cells
    int num_cells = mesh->num_cells();
    Check (absorption.size() == num_cells);

    // return fleck factors
    SP_Fleck fleck(new Fleck_Factors<MT>(mesh));
    Check (fleck->fleck.size() == num_cells);

    // constants needed for calculation of fleck factors
    double T      = 0.0; // temp in keV
    double dedT   = 0.0; // dedT in Jerks/keV
    double volume = 0.0; // volume in cc
    double beta   = 0.0;

    // loop through the cells and assign the fleck factors
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// get material data
	T      = mat_state->get_T(cell);
	dedT   = mat_state->get_dedt(cell);
	volume = mesh->volume(cell);

	Check (T      >=  0.0);
	Check (dedT   >   0.0);
	Check (volume >   0.0);

	// calculate beta (4acT^3/Cv)
	beta = 4.0 * a * T*T*T * volume / dedT;

	// calculate Fleck Factor
	fleck->fleck(cell) = 1.0 / 
	    (1.0 + implicitness * beta * c * delta_t * absorption(cell));
	Check (fleck->fleck(cell) >= 0.0 && fleck->fleck(cell) <= 1.0);
    }   

    Ensure (fleck);
    return fleck;
}

//===========================================================================//
/*!
 * \class Opacity_Builder_Helper<MT,Multigroup_Frequency>
 *
 * \brief Specialization of Opacity_Builder_Helper on Multigroup_Frequency.
 */
//===========================================================================//

template<class MT>
class Opacity_Builder_Helper<MT, Multigroup_Frequency>
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>                                 SP_Mesh;
    typedef rtt_dsxx::SP<Multigroup_Frequency>               SP_Frequency;
    typedef rtt_dsxx::SP<Opacity<MT, Multigroup_Frequency> > SP_Opacity;
    typedef rtt_dsxx::SP<Mat_State<MT> >                     SP_Mat_State;
    typedef rtt_dsxx::SP<Diffusion_Opacity<MT> >             SP_Diff_Opacity;
    typedef std::pair<SP_Opacity, SP_Diff_Opacity>           pair_opacities;
    typedef std::vector<double>                              sf_double;
    typedef typename MT::template CCSF<sf_double>            ccvf;

  public:
    // Build multigroup opacities.
    static pair_opacities build_Opacity(SP_Mesh, SP_Frequency, SP_Mat_State,
					const ccvf &, const ccvf &, 
					const double, const double, bool);
};

//---------------------------------------------------------------------------//
/*!
 * \brief Build Opacity<MT,Multigroup_Opacity> and Diffusion_Opacity<MT>
 * objects. 
 */
template<class MT>
typename Opacity_Builder_Helper<MT,Multigroup_Frequency>::pair_opacities
Opacity_Builder_Helper<MT,Multigroup_Frequency>::build_Opacity(
    SP_Mesh      mesh,
    SP_Frequency frequency,
    SP_Mat_State mat_state,
    const ccvf  &absorption,
    const ccvf  &scattering,
    const double delta_t,
    const double implicitness,
    bool         build_diffusion_opacity)
{
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_cdi::CDI;
    using rtt_dsxx::soft_equiv;
    using rtt_dsxx::SP;
    using std::pair;

    Require (frequency);
    Require (mat_state);
    Require (mesh);
    
    // number of cells and groups
    int num_cells  = mesh->num_cells();
    int num_groups = frequency->get_num_groups();
    Check (num_groups > 0);

    // opacity and diffusion opacity
    SP_Opacity      opacity;
    SP_Diff_Opacity diff_opacity;

    // opacity data
    typename MT::template CCSF<double>    rosseland(mesh);
    typename MT::template CCSF<double>    gray_scat(mesh);
    typename MT::template CCSF<double>    integrated_norm_planck(mesh);
    typename MT::template CCSF<sf_double> emission_group_cdf(mesh);

    // make a Fleck Factors object
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    Check (fleck->fleck.size() == num_cells);

    // variables needed to calculate the fleck factor and integrated Planck
    // and Rosseland functions 
    double planck        = 0.0;
    double dedT          = 0.0;
    double volume        = 0.0;
    double beta          = 0.0;
    double T             = 0.0;
    double b_g           = 0.0;
    double r_g           = 0.0;
    double b_sum         = 0.0;
    double r_sum         = 0.0;
    double sig_p_sum     = 0.0;
    double inv_sig_r_sum = 0.0;
    double inv_sct_r_sum = 0.0;
    double tot_sig_g     = 0.0;

    // group boundaries
    pair<double,double> bounds;

    // loop through cells and build the integrated normalized Planck
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// initialize summation of sigma * b_g over the cell
	sig_p_sum = 0.0;

	// initialize summation of 1/sigma * r_g over the cell
	inv_sig_r_sum = 0.0;
	inv_sct_r_sum = 0.0;

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
	    bounds = frequency->get_group_boundaries(g);

	    // integrate the normalized Planckian and Rosseland functions
	    // over the group
	    CDI::integrate_Rosseland_Planckian_Spectrum(bounds.first, 
							bounds.second,
							mat_state->get_T(cell),
							b_g,
							r_g);

	    // multiply by the absorption opacity and sum
	    sig_p_sum += b_g * absorption(cell)[g-1];

	    // assign to the cdf
	    emission_group_cdf(cell)[g-1] = sig_p_sum;

	    // we calculate the rosseland absorption opacity, not the
	    // rosseland total opacity
	    tot_sig_g = absorption(cell)[g-1];
	    
	    // calculate numerator of Rosseland absorption opacity
	    if (tot_sig_g == 0.0)
	    {
		inv_sig_r_sum += rtt_mc::global::huge;
	    }
	    else
	    {
		inv_sig_r_sum += r_g / tot_sig_g;
	    }

	    // calculate the gray Rosseland scattering opacity
	    tot_sig_g = scattering(cell)[g-1];
	    
	    // calculate numerator of Rosseland scattering opacity
	    if (tot_sig_g == 0.0)
	    {
		inv_sct_r_sum += rtt_mc::global::huge;
	    }
	    else
	    {
		inv_sct_r_sum += r_g / tot_sig_g;
	    }

	}

	// integrate the unnormalized Planckian and Rosseland over the group
	// spectrum 
	CDI::integrate_Rosseland_Planckian_Spectrum(
	    frequency->get_group_boundaries().front(),
	    frequency->get_group_boundaries().back(), 
	    mat_state->get_T(cell),
	    b_sum, r_sum);

	// assign the Planckian integral to the integrated_norm_planck
	integrated_norm_planck(cell) = b_sum;
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

	// build the rosseland absorption opacity
	Check (inv_sig_r_sum > 0.0);
	rosseland(cell) = r_sum / inv_sig_r_sum;

	// build the rosseland scattering opacity
	Check (inv_sct_r_sum > 0.0);
	gray_scat(cell) = r_sum / inv_sct_r_sum;

	// calculate beta (4aT^3/Cv)
	beta = 4.0 * a * T*T*T * volume / dedT;
	
	// calculate Fleck Factor
	fleck->fleck(cell) = 1.0 / 
	    (1.0 + implicitness * beta * c * delta_t * planck); 
	Check (fleck->fleck(cell) >= 0.0 && fleck->fleck(cell) <= 1.0);
    }

    // build the opacity
    opacity = new Opacity<MT, Multigroup_Frequency>(
	frequency, absorption, scattering, fleck, integrated_norm_planck, 
	emission_group_cdf);

    // build the Diffusion_Opacity
    if (build_diffusion_opacity)
    {
	// make the diffusion opacity
	diff_opacity = new Diffusion_Opacity<MT>(fleck, rosseland, gray_scat);
	Ensure (diff_opacity);
	Ensure (diff_opacity->num_cells() == mesh->num_cells());
    }
    Check (opacity);
    Check (opacity->num_cells() == mesh->num_cells());

    // return opacities
    pair_opacities opacities = std::make_pair(opacity, diff_opacity);

    Ensure (opacities.first);
    return opacities;
}

} // end namespace rtt_imc

#endif // rtt_imc_Opacity_Builder_Helper_hh

//---------------------------------------------------------------------------//
//              end of imc/Opacity_Builder_Helper.hh
//---------------------------------------------------------------------------//
