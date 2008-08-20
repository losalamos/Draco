//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Odfmg_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Analytic_Odfmg_Opacity class member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Analytic_Odfmg_Opacity.hh"
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
Analytic_Odfmg_Opacity::Analytic_Odfmg_Opacity(
		const sf_double         &groups,
		const sf_double         &bands,
		const sf_Analytic_Model &models,
		rtt_cdi::Reaction        reaction_in,
		rtt_cdi::Model           model_in)
	: groupBoundaries(groups),
	bandBoundaries(bands),
	group_models(models),
	reaction(reaction_in),
	model(model_in)
{
	Require (reaction == rtt_cdi::TOTAL ||
				reaction == rtt_cdi::ABSORPTION ||
				reaction == rtt_cdi::SCATTERING);
	Require (groupBoundaries.size() - 1 == group_models.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 * 
 * This constructor rebuilds and Analytic_Odfmg_Opacity from a
 * vector<char> that was created by a call to pack().  It can only rebuild
 * Analytic_Model types that have been registered in the
 * rtt_cdi_analytic::Opacity_Models enumeration.
 */
Analytic_Odfmg_Opacity::Analytic_Odfmg_Opacity(
		const sf_char &packed)
{
	// the packed size must be at least 5 integers (number of groups, number of 
	// bands, reaction type, model type, analytic model indicator)
	Require (packed.size() >= 5 * sizeof(int));

	// make an unpacker
	rtt_dsxx::Unpacker unpacker;

	// register the unpacker
	unpacker.set_buffer(packed.size(), &packed[0]);

	// unpack the number of group boundaries
	int ngrp_bounds = 0;
	unpacker >> ngrp_bounds;
	int num_groups  = ngrp_bounds - 1;

	// make the group boundaries and model vectors
	groupBoundaries.resize(ngrp_bounds);
	group_models.resize(num_groups);

	// unpack the group boundaries
	for (int i = 0; i < ngrp_bounds; i++)
		unpacker >> groupBoundaries[i];
	
	// unpack the number of band boundaries
	int nband_bounds = 0;
	unpacker >> nband_bounds;

	// make the group boundaries and model vectors
	bandBoundaries.resize(nband_bounds);

	// unpack the group boundaries
	for (int i = 0; i < nband_bounds; i++)
		unpacker >> bandBoundaries[i];

	// now unpack the models
	std::vector<sf_char> models(num_groups);
	int                  model_size = 0;
	for (int i = 0; i < models.size(); i++)
	{
		// unpack the size of the analytic model
		unpacker >> model_size;
		Check (model_size >= sizeof(int));

		models[i].resize(model_size);

		// unpack the model
		for (int j = 0; j < models[i].size(); j++)
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
	for (int i = 0; i < models.size(); i++)
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

	Ensure (groupBoundaries.size() - 1 == group_models.size());
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
 * (Analytic_Odfmg_Opacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return group opacities (coefficients) in cm^2/g
 *
 */
std::vector< std::vector<double> > Analytic_Odfmg_Opacity::getOpacity( 
		double targetTemperature,
		double targetDensity ) const 
{
	Require (targetTemperature >= 0.0);
	Require (targetDensity >= 0.0);

	const int numBands = getNumBands();
	const int numGroups = getNumGroups();

	// return opacities
	std::vector< std::vector<double> > opacity( numGroups );

	// loop through groups and get opacities
	for (int group = 0; group < opacity.size(); group++)
	{
		Check (group_models[group]);

		opacity[group].resize(numBands);

		// assign the opacity based on the group model to the first band
		opacity[group][0] = group_models[group]->
			calculate_opacity(targetTemperature, targetDensity);

		Check (opacity[group][0] >= 0.0);

		//copy the opacity to the rest of the bands
		for (int band = 1; band < numBands; band++)
		{
			opacity[group][band] = opacity[group][0];
		}
	}

	return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband
 *     opacity 2-D vectors that correspond to the provided vector of
 *     temperatures and a single density value.
 */
std::vector< std::vector< std::vector<double> > >  Analytic_Odfmg_Opacity::getOpacity( 
		const std::vector<double>& targetTemperature,
		double targetDensity ) const
{ 
	std::vector< std::vector< std::vector<double> > > opacity( targetTemperature.size() );

	for ( int i=0; i<targetTemperature.size(); ++i )
	{
		opacity[i] = getOpacity(targetTemperature[i], targetDensity);
	}
	return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband
 *     opacity 2-D vectors that correspond to the provided
 *     temperature and a vector of density values.
 */
std::vector< std::vector< std::vector<double> > >  Analytic_Odfmg_Opacity::getOpacity( 
		double targetTemperature,
		const std::vector<double>& targetDensity ) const
{ 
	std::vector< std::vector< std::vector<double> > > opacity( targetDensity.size() );

	//call our regular getOpacity function for every target density
	for ( int i=0; i<targetDensity.size(); ++i )
	{
		opacity[i] = getOpacity(targetTemperature, targetDensity[i]);
	}
	return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic odfmg opacity.
 *
 * This function will pack up the Analytic_Mulitgroup_Opacity into a char
 * array (represented by a vector<char>).  The Analytic_Opacity_Model derived
 * class must have a pack function; this is enforced by the virtual
 * Analytic_Opacity_Model base class.
 */
Analytic_Odfmg_Opacity::sf_char Analytic_Odfmg_Opacity::pack() const
{
	Require (groupBoundaries.size() - 1 == group_models.size());

	// make a packer
	rtt_dsxx::Packer packer;

	// first pack up models
	std::vector<sf_char> models(group_models.size());
	int                  num_bytes_models = 0;

	// loop through and pack up the models
	for (int i = 0; i < models.size(); i++)
	{
		Check (group_models[i]);

		models[i]         = group_models[i]->pack();
		num_bytes_models += models[i].size();
	}

	// now add up the total size; number of groups + 1 int for number of
	// groups, number of bands + 1 int for number of
	// bands, number of models + size in each model + models, 1 int for
	// reaction type, 1 int for model type
	int size = 4 * sizeof(int) + models.size() * sizeof(int) + 
		groupBoundaries.size() * sizeof(double) +
	   bandBoundaries.size()  * sizeof(double) + num_bytes_models;

	// make a char array
	sf_char packed(size);

	// set the buffer
	packer.set_buffer(size, &packed[0]);

	// pack the number of groups and group boundaries
	packer << static_cast<int>(groupBoundaries.size());
	for (int i = 0; i < groupBoundaries.size(); i++)
		packer << groupBoundaries[i];

	// pack the number of bands and band boundaries
	packer << static_cast<int>(bandBoundaries.size());
	for (int i = 0; i < bandBoundaries.size(); i++)
		packer << bandBoundaries[i];

	// pack each models size and data
	for (int i = 0; i < models.size(); i++)
	{
		// pack the size of this model
		packer << static_cast<int>(models[i].size());

		// now pack the model data
		for (int j = 0; j < models[i].size(); j++)
			packer << models[i][j];
	}

	// now pack the reaction and model type
	packer << static_cast<int>(reaction) << static_cast<int>(model);

	Ensure (packer.get_ptr() == &packed[0] + size);
	return packed;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Analytic_Odfmg_Opacity.cc
//---------------------------------------------------------------------------//
