//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Gray_Opacity.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 24 13:13:46 2001
 * \brief  Analytic_Gray_Opacity member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Analytic_Gray_Opacity.hh"
#include <cmath>

namespace rtt_cdi_analytic
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic opacity model.
 *
 * This constructor builds an opacity model defined by the
 * rtt_cdi_analytic::Analytic_Opacity_Model derived class argument.
 *
 * The reaction type for this instance of the class is determined by the
 * rtt_cdi::Reaction argument.
 *
 * \param analytic_model_in rtt_dsxx::SP to a derived
 * rtt_cdi_analytic::Analytic_Opacity_Model object
 *
 * \param reaction_in rtt_cdi::Reaction type (enumeration)
 *
 */
Analytic_Gray_Opacity::Analytic_Gray_Opacity(SP_Analytic_Model model_in,
					     rtt_cdi::Reaction reaction_in)
    : analytic_model(model_in),
      reaction(reaction_in)
{
    Require (reaction == rtt_cdi::TOTAL ||
	     reaction == rtt_cdi::ABSORPTION ||
	     reaction == rtt_cdi::SCATTERING);

    Ensure (analytic_model);
}

//---------------------------------------------------------------------------//
// OPACITY INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return a scalar opacity given a scalar temperature and density. 
 *
 * Given a scalar temperature and density, return an opacity for the reaction
 * type specified by the constructor.  The analytic opacity function is
 * specified in the constructor (Analytic_Gray_Opacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return opacity (coefficient) in cm^2/g
 *
 */
double Analytic_Gray_Opacity::getOpacity(double temperature,
					 double density) const
{
    Require (temperature >= 0.0);
    Require (density >= 0.0);

    // define return opacity
    double opacity = analytic_model->calculate_opacity(temperature, density);

    Ensure (opacity >= 0.0);
    return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a field of opacities given a field of temperatures and a
 * scalar density. 
 *
 * Given a field of temperatures and a scalar density, return an opacity
 * field for the reaction type specified by the constructor.  The analytic
 * opacity function is specified in the constructor
 * (Analytic_Gray_Opacity()).  The returned opacity field has the same number
 * of elements as the temperature field.
 *
 * The field type for temperatures is an std::vector.
 *
 * \param temperature std::vector of material temperatures in keV
 * \param density material density in g/cm^3
 * \return std::vector of opacities (coefficients) in cm^2/g
 */
Analytic_Gray_Opacity::sf_double
Analytic_Gray_Opacity::getOpacity(const sf_double &temperature,
				  double density) const
{
    Require (density >= 0.0);

    // define the return opacity field (same size as temperature field)
    sf_double opacity(temperature.size(), 0.0);

    // define an opacity iterator
    sf_double::iterator sig = opacity.begin();

    // loop through temperatures and solve for opacity
    for (sf_double::const_iterator T = temperature.begin(); 
	 T != temperature.end(); T++, sig++)
    {
	Check (*T >= 0.0);

	// define opacity
	*sig = analytic_model->calculate_opacity(*T, density);

	Check (*sig >= 0.0);
    }

    return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a field of opacities given a field of densities and a scalar
 * temperature.
 *
 * Given a field of densities and a scalar temperature, return an opacity
 * field for the reaction type specified by the constructor.  The analytic
 * opacity function is specified in the constructor
 * (Analytic_Gray_Opacity()).  The returned opacity field has the same number
 * of elements as the density field.
 *
 * The field type for densities is an std::vector.
 *
 * \param temperature material temperature in keV
 * \param density std::vector of material densities in g/cc
 * \return std::vector of opacities (coefficients) in cm^2/g
 */
Analytic_Gray_Opacity::sf_double
Analytic_Gray_Opacity::getOpacity(double temperature,
				  const sf_double &density) const
{
    Require (temperature >= 0.0);

    // define the return opacity field (same size as density field)
    sf_double opacity(density.size(), 0.0);

    // define an opacity iterator
    sf_double::iterator sig = opacity.begin();

    // loop through densities and solve for opacity
    for (sf_double::const_iterator rho = density.begin();
	 rho != density.end(); rho++, sig++)
    {
	Check (*rho >= 0.0);

	// define opacity
	*sig = analytic_model->calculate_opacity(temperature, *rho);

	Check (*sig >= 0.0);
    }

    return opacity;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Analytic_Gray_Opacity.cc
//---------------------------------------------------------------------------//
