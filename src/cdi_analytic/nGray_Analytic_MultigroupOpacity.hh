//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/nGray_Analytic_MultigroupOpacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  nGray_Analytic_MultigroupOpacity class definition.
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_nGray_Analytic_MultigroupOpacity_hh__
#define __cdi_analytic_nGray_Analytic_MultigroupOpacity_hh__

#include "Analytic_Models.hh"
#include "cdi/MultigroupOpacity.hh"
#include "cdi/OpacityCommon.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_cdi_analytic
{

//===========================================================================//
/*!
 * \class nGray_Analytic_MultigroupOpacity
 *
 * \brief Derived rtt_cdi::MultigroupOpacity class for analytic opacities.
 *
 * The nGray_Analytic_MultigroupOpacity class is a derived
 * rtt_cdi::MultigroupOpacity class.  It provides analytic opacity data. The
 * specific analytic opacity model is derived from the
 * rtt_cdi_analytic::Analytic_Opacity_Model base class.  Several pre-built
 * derived classes are provided in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it
 * conforms to the rtt_cdi_analytic::Analytic_Opacity_Model interface.  This
 * interface consists of a single function,
 * Analytic_Opacity_Model::calculate_opacity().
 *
 * An rtt_cdi_analytic::Analytic_Opacity_Model is provided for each analytic
 * energy group.  In other words, the user specifies the group structure for
 * the multigroup data and applies a Analytic_Model for each group.  This way
 * a different rtt_cdi_analytic::Analytic_Model can be used in each group.
 * For example, the client could choose to use a
 * rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model for the low frequency
 * groups and a rtt_cdi_analytic::Constant_Analytic_Opacity_Model for the
 * high frequency groups.
 *
 * Note that opacities are returned in units of cm^2/g. Thus the resultant
 * opacity must be multiplied by density to get units of 1/cm.  See the
 * documentation in rtt_cdi_analytic::Analytic_Model for more info.
 *
 * The constructors take a rtt_cdi::Reaction argument to determine the
 * reaction type.  The enumeration rtt_cdi::Reaction can have the value
 * TOTAL, ABSORPTION, or SCATTERING.
 *
 * The default rtt_cdi::Model for an Analytic_Gray_Opacity is
 * rtt_cdi::ANALYTIC.  However, this can be overridden in the constructor.
 *
 * This class conforms to the interface specified by
 * rtt_cdi::MultigroupOpacity and can be used with rtt_cdi::CDI to get
 * analytic opacities.
 *
 * \example cdi_analytic/test/tstnGray_Analytic_MultigroupOpacity.cc
 *
 * Example usage of nGray_Analytic_MultigroupOpacity, Analytic_Opacity_Model, and
 * their incorporation into rtt_cdi::CDI.
 */ 
// revision history:
// -----------------
// 0) original
// 1) 06-MAR-03 : added model specification; default is still ANALYTIC
// 
//===========================================================================//

 class nGray_Analytic_MultigroupOpacity : public rtt_cdi::MultigroupOpacity
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Analytic_Opacity_Model>       SP_Analytic_Model;
    typedef rtt_dsxx::SP<const Analytic_Opacity_Model> const_Model; 
    typedef std::vector<SP_Analytic_Model>             sf_Analytic_Model;
    typedef std::vector<double>                        sf_double;
    typedef std::vector<sf_double>                     vf_double;
    typedef std::string                                std_string;
    typedef std::vector<char>                          sf_char;
    
  private:
    // Group structure.
    sf_double group_boundaries;

    // Analytic models for each group.
    sf_Analytic_Model group_models;

    // Reaction model.
    rtt_cdi::Reaction reaction;

    // CDI model.
    rtt_cdi::Model model;

  public:
    // Constructor.
    nGray_Analytic_MultigroupOpacity(const sf_double &, 
				const sf_Analytic_Model &,
				rtt_cdi::Reaction,
				rtt_cdi::Model = rtt_cdi::ANALYTIC);

    // Constructor for packed Analytic_Multigroup_Opacities
    explicit nGray_Analytic_MultigroupOpacity(const sf_char &);

    // >>> ACCESSORS
    const_Model get_Analytic_Model(size_t g) const { return group_models[g-1]; }

    // >>> INTERFACE SPECIFIED BY rtt_cdi::MultigroupOpacity

    // Get the group opacities.
    sf_double getOpacity(double, double) const;

    // Get the group opacity fields given a field of temperatures.
    vf_double getOpacity(const sf_double &, double) const;

    // Get the group opacity fields given a field of densities.
    vf_double getOpacity(double, const sf_double &) const;

    //! Query to see if data is in tabular or functional form (false).
    bool data_in_tabular_form() const { return false; }

    //! Query to get the reaction type.
    rtt_cdi::Reaction getReactionType() const { return reaction; }

    //! Query for model type.
    rtt_cdi::Model getModelType() const { return model; }

    // Return the energy policy (gray).
    inline std_string getEnergyPolicyDescriptor() const;

    // Get the data description of the opacity.
    inline std_string getDataDescriptor() const;

    // Get the name of the associated data file.
    inline std_string getDataFilename() const;

    //! Get the temperature grid (size 0 for function-based analytic data).
    sf_double getTemperatureGrid() const { return sf_double(0); }

    //! Get the density grid (size 0 for function-based analytic data).
    sf_double getDensityGrid() const { return sf_double(0); }

    //! Get the group boundaries (keV) of the multigroup set.
    sf_double getGroupBoundaries() const { return group_boundaries; }

    //! Get the size of the temperature grid (size 0).
    size_t getNumTemperatures() const { return 0; }

    //! Get the size of the density grid (size 0).
    size_t getNumDensities() const { return 0; }

    //! Get the number of frequency group boundaries.
    size_t getNumGroupBoundaries() const { return group_boundaries.size(); }

    //! Get the number of frequency group boundaries.
    size_t getNumGroups() const { return group_boundaries.size() - 1; }

    // Pack the nGray_Analytic_MultigroupOpacity into a character string.
    sf_char pack() const;

	/*!
	 * \brief Returns the general opacity model type, defined in OpacityCommon.hh
	 *
	 * Since this is an analytic model, return 1 (rtt_cdi::ANALYTIC_TYPE)
	 */
	rtt_cdi::OpacityModelType getOpacityModelType() const {
		return rtt_cdi::ANALYTIC_TYPE;
	}
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the energy policy descriptor (mg). 
 */
nGray_Analytic_MultigroupOpacity::std_string 
nGray_Analytic_MultigroupOpacity::getEnergyPolicyDescriptor() const 
{
    return std_string("mg");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a string describing the opacity model.
 */
nGray_Analytic_MultigroupOpacity::std_string
nGray_Analytic_MultigroupOpacity::getDataDescriptor() const
{
    std_string descriptor;

    if (reaction == rtt_cdi::TOTAL)
	descriptor = "Analytic Multigroup Total";
    else if (reaction == rtt_cdi::ABSORPTION)
	descriptor = "Analytic Multigroup Absorption";
    else if (reaction == rtt_cdi::SCATTERING)
	descriptor = "Analytic Multigroup Scattering";
    else
    {
	Insist (0, "Invalid analytic multigroup model opacity!");
    }

    return descriptor;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return NULL string for the data filename.
 */
nGray_Analytic_MultigroupOpacity::std_string
nGray_Analytic_MultigroupOpacity::getDataFilename() const 
{
    return std_string();
}

} // end namespace rtt_cdi_analytic

#endif              // __cdi_analytic_nGray_Analytic_MultigroupOpacity_hh__

//---------------------------------------------------------------------------//
//            end of cdi_analytic/nGray_Analytic_MultigroupOpacity.hh
//---------------------------------------------------------------------------//
