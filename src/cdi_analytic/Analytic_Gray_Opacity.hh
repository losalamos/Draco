//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Gray_Opacity.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 24 13:13:46 2001
 * \brief  Analytic_Gray_Opacity class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Analytic_Gray_Opacity_hh__
#define __cdi_analytic_Analytic_Gray_Opacity_hh__

#include "Analytic_Models.hh"
#include "cdi/GrayOpacity.hh"
#include "cdi/OpacityCommon.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_cdi_analytic
{
 
//===========================================================================//
/*!
 * \class Analytic_Gray_Opacity
 *
 * \brief Derived rtt_cdi::GrayOpacity class for analytic opacities.
 *
 * The Analytic_Gray_Opacity class is a derived rtt_cdi::GrayOpacity class.
 * It provides analytic opacity data. The specific analytic opacity model is
 * derived from the rtt_cdi_analytic::Analytic_Opacity_Model base class.
 * Several pre-built derived classes are provided in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it
 * conforms to the rtt_cdi_analytic::Analytic_Opacity_Model interface.  This
 * interface consists of a single function,
 * Analytic_Opacity_Model::calculate_opacity().
 *
 * Note that opacities are returned in units of cm^2/g. Thus the resultant
 * opacity must be multiplied by density to get units of 1/cm.  See the
 * documentation in rtt_cdi_analytic::Analytic_Model for more info.
 *
 * The constructors take a rtt_cdi::Reaction argument to determine the
 * reaction type.  The enumeration rtt_cdi::Reaction can have the value
 * TOTAL, ABSORPTION, or SCATTERING.
 *
 * This class conforms to the interface specified by rtt_cdi::GrayOpacity and
 * can be used with rtt_cdi::CDI to get analytic opacities.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Analytic_Gray_Opacity : public rtt_cdi::GrayOpacity
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Analytic_Opacity_Model> SP_Analytic_Model;
    typedef std::vector<double>                  sf_double;
    typedef std::string                          std_string;
    
  private:
    // Analytic opacity model.
    SP_Analytic_Model analytic_model;

    // Reaction model.
    rtt_cdi::Reaction reaction;
    
  public:
    // Constructor.
    Analytic_Gray_Opacity(SP_Analytic_Model, rtt_cdi::Reaction);

    // >>> INTERFACE SPECIFIED BY rtt_cdi::GrayOpacity

    // Get an opacity.
    double getOpacity(double, double) const;

    // Get an opacity field given a field of temperatures.
    sf_double getOpacity(const sf_double &, double) const;

    // Get an opacity field given a field of densities.
    sf_double getOpacity(double, const sf_double &) const;

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

    //! Get the size of the temperature grid (size 0).
    int getNumTemperatures() const { return 0; }

    //! Get the size of the density grid (size 0).
    int getNumDensities() const { return 0; }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the energy policy descriptor (gray for
 * Analytic_Gray_Opacity). 
 */
Analytic_Gray_Opacity::std_string 
Analytic_Gray_Opacity::getEnergyPolicyDescriptor() const 
{
    return std_string("gray");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a string describing the opacity model.
 */
Analytic_Gray_Opacity::std_string
Analytic_Gray_Opacity::getDataDescriptor() const
{
    std_string descriptor;

    if (reaction == rtt_cdi::TOTAL)
	descriptor = "Analytic Gray Total";
    else if (reaction == rtt_cdi::ABSORPTION)
	descriptor = "Analytic Gray Absorption";
    else if (reaction == rtt_cdi::SCATTERING)
	descriptor = "Analytic Gray Scattering";
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
Analytic_Gray_Opacity::std_string
Analytic_Gray_Opacity::getDataFilename() const 
{
    return std_string(0);
}

} // end namespace rtt_cdi_analytic

#endif                          // __cdi_analytic_Analytic_Gray_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_analytic/Analytic_Gray_Opacity.hh
//---------------------------------------------------------------------------//
