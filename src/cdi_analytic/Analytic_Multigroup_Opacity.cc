//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Multigroup_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Analytic_Multigroup_Opacity class member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Analytic_Multigroup_Opacity.hh"

namespace rtt_cdi_analytic
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic multigroup opacity model.
 *
 * This constructor builds an opacity model defined by the
 * rtt_cdi_analytic::Analytic_Opacity_Model derived class argument.
 *
 * The reaction type for this instance of the class is determined by the
 * rtt_cdi::Reaction argument.
 *
 * The group structure (in keV) must be provided by the groups argument.  The
 * number of Analytic_Opacity_Model objects given in the models argument must
 * be equal to the number of groups.
 *
 * \param groups vector containing the group boundaries in keV from lowest to
 * highest
 *
 * \param models vector containing SPs to Analytic_Model derived types for
 * each group, the size should be groups.size() - 1
 *
 * \param reaction_in rtt_cdi::Reaction type (enumeration)
 *
 */
Analytic_Multigroup_Opacity::Analytic_Multigroup_Opacity(
    const sf_double &groups,
    const sf_Analytic_Model &models,
    rtt_cdi::Reaction reaction_in)
    : group_boundaries(groups),
      group_models(models),
      reaction(reaction_in)
{
    Require (reaction == rtt_cdi::TOTAL ||
	     reaction == rtt_cdi::ABSORPTION ||
	     reaction == rtt_cdi::SCATTERING);
    Require (group_boundaries.size() - 1 == group_models.size());
}

//---------------------------------------------------------------------------//
// OPACITY INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the group opacities given a scalar temperature and density. 
 *
 * Given a scalar temperature and density, return the group opacities
 * (vector<double>) for the reaction type specified by the constructor.  The
 * analytic opacity model is specified in the constructor
 * (Analytic_Multigroup_Opacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return group opacities (coefficients) in cm^2/g
 *
 */
Analytic_Multigroup_Opacity::sf_double
Analytic_Multigroup_Opacity::getOpacity(double temperature, 
					double density) const
{
    Require (temperature >= 0.0);
    Require (density >= 0.0);

    // return opacities
    sf_double opacities(group_models.size(), 0.0);

    // loop through groups and get opacities
    for (int i = 0; i < opacities.size(); i++)
    {
	Check (group_models[i]);

	// assign the opacity based on the group model
	opacities[i] = group_models[i]->
	    calculate_opacity(temperature, density);

	Check (opacities[i] >= 0.0);
    }
    
    return opacities;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Return a vector of multigroup opacities given a vector of
 * temperatures and a scalar density.
 *
 * Given a field of temperatures and a scalar density, return a vector of
 * multigroup opacities (vector<vector<double>>) for the reaction type
 * specified by the constructor.  The analytic opacity model is specified in
 * the constructor (Analytic_Gray_Opacity()).  The returned opacity field is
 * indexed [num_temperatures][num_groups].
 *
 * The field type for temperatures is an std::vector.
 *
 * \param temperature std::vector of material temperatures in keV 
 *
 * \param density material density in g/cm^3
 *
 * \return std::vector<std::vector> of multigroup opacities (coefficients) in
 * cm^2/g indexed [temperature][group]
 */
Analytic_Multigroup_Opacity::vf_double
Analytic_Multigroup_Opacity::getOpacity(const sf_double &temperature,
					double density) const
{
    Require (density >= 0.0);

    // define the return opacity field (same size as temperature field); each
    // entry has a vector sized by the number of groups
    vf_double opacities(temperature.size(), 
			sf_double(group_models.size(), 0.0));

    // loop through temperatures and solve for opacity
    for (int i = 0; i < opacities.size(); i++)
    {
	Check (temperature[i] >= 0.0);

	// loop through groups
	for (int j = 0; j < opacities[i].size(); j++)
	{
	    Check (group_models[j]);

	    // assign the opacity based on the group model
	    opacities[i][j] = group_models[j]->
		calculate_opacity(temperature[i], density);

	    Check (opacities[i][j] >= 0.0);
	}
    }

    return opacities;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Return a vector of multigroup opacities given a vector of
 * density and a scalar temperature.
 *
 * Given a field of densities and a scalar temperature, return a vector of
 * multigroup opacities (vector<vector<double>>) for the reaction type
 * specified by the constructor.  The analytic opacity model is specified in
 * the constructor (Analytic_Gray_Opacity()).  The returned opacity field is
 * indexed [num_density][num_groups].
 *
 * The field type for density is an std::vector.
 *
 * \param temperature in keV 
 *
 * \param density std::vector of material densities in g/cm^3
 *
 * \return std::vector<std::vector> of multigroup opacities (coefficients) in
 * cm^2/g indexed [density][group]
 */
Analytic_Multigroup_Opacity::vf_double
Analytic_Multigroup_Opacity::getOpacity(double temperature,
					const sf_double &density) const
{
    Require (temperature >= 0.0);

    // define the return opacity field (same size as density field); each
    // entry has a vector sized by the number of groups
    vf_double opacities(density.size(), 
			sf_double(group_models.size(), 0.0));

    // loop through densities and solve for opacity
    for (int i = 0; i < opacities.size(); i++)
    {
	Check (density[i] >= 0.0);

	// loop through groups
	for (int j = 0; j < opacities[i].size(); j++)
	{
	    Check (group_models[j]);

	    // assign the opacity based on the group model
	    opacities[i][j] = group_models[j]->
		calculate_opacity(temperature, density[i]);

	    Check (opacities[i][j] >= 0.0);
	}
    }

    return opacities;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Analytic_Multigroup_Opacity.cc
//---------------------------------------------------------------------------//
