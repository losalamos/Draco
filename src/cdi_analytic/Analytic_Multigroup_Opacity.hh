//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Multigroup_Opacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Analytic_Multigroup_Opacity class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Analytic_Multigroup_Opacity_hh__
#define __cdi_analytic_Analytic_Multigroup_Opacity_hh__

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
 * \class Analytic_Multigroup_Opacity
 *
 * \brief Derived rtt_cdi::MultigroupOpacity class for analytic opacities.
 *
 * The Analytic_Multigroup_Opacity class is a derived
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
 * This class conforms to the interface specified by
 * rtt_cdi::MultigroupOpacity and can be used with rtt_cdi::CDI to get
 * analytic opacities.
 *
 */
/*!
 *
 * \example cdi_analytic/test/tstAnalytic_Multigroup_Opacity.cc
 *
 * Example usage of Analytic_Multigroup_Opacity, Analytic_Opacity_Model, and
 * their incorporation into rtt_cdi::CDI.
 *
 */ 
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

 class Analytic_Multigroup_Opacity : public rtt_cdi::MultigroupOpacity
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Analytic_Opacity_Model> SP_Analytic_Model;
    typedef std::vector<SP_Analytic_Model>       sf_Analytic_Model;
    typedef std::vector<double>                  sf_double;
    typedef std::vector<sf_double>               vf_double;
    typedef std::string                          std_string;
    typedef std::vector<char>                    sf_char;
    
  private:
    // Group structure.
    sf_double group_boundaries;

    // Analytic models for each group.
    sf_Analytic_Model group_models;

    // Reaction model.
    rtt_cdi::Reaction reaction;

  public:
    // Constructor.
    Analytic_Multigroup_Opacity(const sf_double &, 
				const sf_Analytic_Model &,
				rtt_cdi::Reaction);

    // Constructor for packed Analytic_Multigroup_Opacities
    Analytic_Multigroup_Opacity(const sf_char &);

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
    rtt_cdi::Model getModelType() const { return rtt_cdi::ANALYTIC; }

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
    int getNumTemperatures() const { return 0; }

    //! Get the size of the density grid (size 0).
    int getNumDensities() const { return 0; }

    //! Get the number of frequency group boundaries.
    int getNumGroupBoundaries() const { return group_boundaries.size(); }

    //! Get the number of frequency group boundaries.
    int getNumGroups() const { return group_boundaries.size() - 1; }

    // Pack the Analytic_Multigroup_Opacity into a character string.
    sf_char pack() const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the energy policy descriptor (mg). 
 */
Analytic_Multigroup_Opacity::std_string 
Analytic_Multigroup_Opacity::getEnergyPolicyDescriptor() const 
{
    return std_string("mg");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a string describing the opacity model.
 */
Analytic_Multigroup_Opacity::std_string
Analytic_Multigroup_Opacity::getDataDescriptor() const
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
	Insist (0, "Invalid analytic gray model opacity!");
    }

    return descriptor;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return NULL string for the data filename.
 */
Analytic_Multigroup_Opacity::std_string
Analytic_Multigroup_Opacity::getDataFilename() const 
{
    return std_string();
}

} // end namespace rtt_cdi_analytic

#endif              // __cdi_analytic_Analytic_Multigroup_Opacity_hh__

//---------------------------------------------------------------------------//
//                   end of cdi_analytic/Analytic_Multigroup_Opacity.hh
//---------------------------------------------------------------------------//
