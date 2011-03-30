//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/nGray_Analytic_MultigroupOpacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  nGray_Analytic_MultigroupOpacity class member definitions.
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "nGray_Analytic_MultigroupOpacity.hh"
#include "ds++/Packing_Utils.hh"

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
nGray_Analytic_MultigroupOpacity::nGray_Analytic_MultigroupOpacity(
    const sf_double         &groups,
    const sf_Analytic_Model &models,
    rtt_cdi::Reaction        reaction_in,
    rtt_cdi::Model           model_in)
    : group_boundaries(groups),
      group_models(models),
      reaction(reaction_in),
      model(model_in)
{
    Require (reaction == rtt_cdi::TOTAL ||
	     reaction == rtt_cdi::ABSORPTION ||
	     reaction == rtt_cdi::SCATTERING);
    Require (group_boundaries.size() - 1 == group_models.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 * 
 * This constructor rebuilds and nGray_Analytic_MultigroupOpacity from a
 * vector<char> that was created by a call to pack().  It can only rebuild
 * Analytic_Model types that have been registered in the
 * rtt_cdi_analytic::Opacity_Models enumeration.
 */
nGray_Analytic_MultigroupOpacity::nGray_Analytic_MultigroupOpacity(
    const sf_char &packed)
    : group_boundaries(),
      group_models(),
      reaction(),
      model()
{
    // the packed size must be at least 4 integers (number of groups,
    // reaction type, model type, analytic model indicator)
    Require (packed.size() >= 4 * sizeof(int));

    // make an unpacker
    rtt_dsxx::Unpacker unpacker;
    
    // register the unpacker
    unpacker.set_buffer(packed.size(), &packed[0]);

    // unpack the number of group boundaries
    int ngrp_bounds = 0;
    unpacker >> ngrp_bounds;
    int num_groups  = ngrp_bounds - 1;
    
    // make the group boundaries and model vectors
    group_boundaries.resize(ngrp_bounds);
    group_models.resize(num_groups);

    // unpack the group boundaries
    for (int i = 0; i < ngrp_bounds; i++)
	unpacker >> group_boundaries[i];

    // now unpack the models
    std::vector<sf_char> models(num_groups);
    int                  model_size = 0;
    for (size_t i = 0; i < models.size(); ++i)
    {
	// unpack the size of the analytic model
	unpacker >> model_size;
	Check (static_cast<size_t>(model_size) >= sizeof(int));

	models[i].resize(model_size);
	
	// unpack the model
	for (size_t j = 0; j < models[i].size(); ++j)
	    unpacker >> models[i][j];
    }

    // unpack the reaction and model type
    int reaction_int, model_int;
    unpacker >> reaction_int >> model_int;
    Check (unpacker.get_ptr() == &packed[0] + packed.size());

    // assign the reaction and model type
    reaction = static_cast<rtt_cdi::Reaction>(reaction_int);
    model    = static_cast<rtt_cdi::Model>(model_int);
    Check (reaction == rtt_cdi::ABSORPTION ||
	   reaction == rtt_cdi::SCATTERING ||
	   reaction == rtt_cdi::TOTAL);

    // now rebuild the analytic models
    int indicator = 0;
    for (size_t i = 0; i < models.size(); ++i)
    {
	// reset the buffer
	unpacker.set_buffer(models[i].size(), &models[i][0]);

	// get the indicator for this model (first packed datum)
	unpacker >> indicator;

	// now determine which analytic model we need to build
	if (indicator == CONSTANT_ANALYTIC_OPACITY_MODEL)
	{
	    group_models[i] = new Constant_Analytic_Opacity_Model(
		models[i]);
	}
	else if (indicator == POLYNOMIAL_ANALYTIC_OPACITY_MODEL)
	{
	    group_models[i] = new Polynomial_Analytic_Opacity_Model(
		models[i]);
	}
	else
	{
	    Insist (0, "Unregistered analytic opacity model!");
	}

	Ensure (group_models[i]);
    }
    
    Ensure (group_boundaries.size() - 1 == group_models.size());
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
 * (nGray_Analytic_MultigroupOpacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return group opacities (coefficients) in cm^2/g
 *
 */
nGray_Analytic_MultigroupOpacity::sf_double
nGray_Analytic_MultigroupOpacity::getOpacity(double temperature, 
					double density) const
{
    Require (temperature >= 0.0);
    Require (density >= 0.0);

    // return opacities
    sf_double opacities(group_models.size(), 0.0);

    // loop through groups and get opacities
    for (size_t i = 0; i < opacities.size(); ++i)
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
nGray_Analytic_MultigroupOpacity::vf_double
nGray_Analytic_MultigroupOpacity::getOpacity(const sf_double &temperature,
					double density) const
{
    Require (density >= 0.0);

    // define the return opacity field (same size as temperature field); each
    // entry has a vector sized by the number of groups
    vf_double opacities(temperature.size(), 
			sf_double(group_models.size(), 0.0));

    // loop through temperatures and solve for opacity
    for (size_t i = 0; i < opacities.size(); ++i)
    {
	Check (temperature[i] >= 0.0);

	// loop through groups
	for (size_t j = 0; j < opacities[i].size(); ++j)
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
nGray_Analytic_MultigroupOpacity::vf_double
nGray_Analytic_MultigroupOpacity::getOpacity(double temperature,
					const sf_double &density) const
{
    Require (temperature >= 0.0);

    // define the return opacity field (same size as density field); each
    // entry has a vector sized by the number of groups
    vf_double opacities(density.size(), 
			sf_double(group_models.size(), 0.0));

    // loop through densities and solve for opacity
    for (size_t i = 0; i < opacities.size(); ++i)
    {
	Check (density[i] >= 0.0);

	// loop through groups
	for (size_t j = 0; j < opacities[i].size(); ++j)
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

//---------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic multigroup opacity.
 *
 * This function will pack up the Analytic_Mulitgroup_Opacity into a char
 * array (represented by a vector<char>).  The Analytic_Opacity_Model derived
 * class must have a pack function; this is enforced by the virtual
 * Analytic_Opacity_Model base class.
 */
nGray_Analytic_MultigroupOpacity::sf_char nGray_Analytic_MultigroupOpacity::pack() const
{
    Require (group_boundaries.size() - 1 == group_models.size());

    // make a packer
    rtt_dsxx::Packer packer;

    // first pack up models
    std::vector<sf_char> models(group_models.size());
    int                  num_bytes_models = 0;

    // loop through and pack up the models
    for (size_t i = 0; i < models.size(); ++i)
    {
	Check (group_models[i]);

	models[i]         = group_models[i]->pack();
	num_bytes_models += models[i].size();
    }

    // now add up the total size; number of groups + 1 int for number of
    // groups, number of models + size in each model + models, 1 int for
    // reaction type, 1 int for model type
    int size = (3 + models.size()) * sizeof(int) + 
	group_boundaries.size() * sizeof(double) + num_bytes_models;

    // make a char array
    sf_char packed(size);

    // set the buffer
    packer.set_buffer(size, &packed[0]);

    // pack the number of groups and group boundaries
    packer << static_cast<int>(group_boundaries.size());
    for (size_t i = 0; i < group_boundaries.size(); ++i)
	packer << group_boundaries[i];

    // pack each models size and data
    for (size_t i = 0; i < models.size(); ++i)
    {
	// pack the size of this model
	packer << static_cast<int>(models[i].size());
	
	// now pack the model data
	for (size_t j = 0; j < models[i].size(); ++j)
	    packer << models[i][j];
    }

    // now pack the reaction and model type
    packer << static_cast<int>(reaction) << static_cast<int>(model);

    Ensure (packer.get_ptr() == &packed[0] + size);
    return packed;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of nGray_Analytic_MultigroupOpacity.cc
//---------------------------------------------------------------------------//
