//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_MultigroupOpacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Analytic_MultigroupOpacity class definition.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Analytic_MultigroupOpacity_hh__
#define __cdi_analytic_Analytic_MultigroupOpacity_hh__

#include "Analytic_Models.hh"
#include "cdi/MultigroupOpacity.hh"
#include <memory>

namespace rtt_cdi_analytic {

//===========================================================================//
/*!
 * \class Analytic_MultigroupOpacity
 *
 * \brief Derived rtt_cdi::MultigroupOpacity class for analytic opacities.
 *
 * The Analytic_MultigroupOpacity class is an abstract class derived from
 * rtt_cdi::MultigroupOpacity class.  It provides analytic opacity data. The
 * specific analytic opacity model is specified by concrete classes derived
 * from Analytic_MultigroupOpacity.
 *
 * Note that opacities are returned in units of cm^2/g. Thus the resultant
 * opacity must be multiplied by density to get units of 1/cm.  See the
 * documentation in rtt_cdi_analytic::Analytic_Model for more info.
 *
 * The constructors take a rtt_cdi::Reaction argument to determine the
 * reaction type.  The enumeration rtt_cdi::Reaction can have the value
 * TOTAL, ABSORPTION, or SCATTERING.
 *
 * The default rtt_cdi::Model for an Analytic_MultigroupOpacity is
 * rtt_cdi::ANALYTIC.  However, this can be overridden in the constructor.
 *
 * \sa cdi_analytic/nGray_Analytic_MultigroupOpacity.hh
 * Example usage of Analytic_MultigroupOpacity.
 */
//===========================================================================//

class Analytic_MultigroupOpacity : public rtt_cdi::MultigroupOpacity {
public:
  // Useful typedefs.
  typedef std::vector<double> sf_double;
  typedef std::vector<sf_double> vf_double;
  typedef std::string std_string;
  typedef std::vector<char> sf_char;

private:
  // Group structure.
  sf_double group_boundaries;

  // Reaction model.
  rtt_cdi::Reaction reaction;

  // CDI model.
  rtt_cdi::Model model;

protected:
  // Constructor.
  Analytic_MultigroupOpacity(const sf_double &groups,
                             rtt_cdi::Reaction reaction_in,
                             rtt_cdi::Model model_in = rtt_cdi::ANALYTIC);

  // Constructor for packed Analytic_Multigroup_Opacities
  explicit Analytic_MultigroupOpacity(const sf_char &packed);

  // Get the packed size of the object
  unsigned packed_size() const;

public:
  // >>> ACCESSORS

  // >>> INTERFACE SPECIFIED BY rtt_cdi::MultigroupOpacity

  virtual ~Analytic_MultigroupOpacity() { /*empty*/
  }

  // Get the group opacities.
  virtual sf_double getOpacity(double, double) const = 0;

  // Get the group opacity fields given a field of temperatures.
  virtual vf_double getOpacity(const sf_double &, double) const = 0;

  // Get the group opacity fields given a field of densities.
  virtual vf_double getOpacity(double, const sf_double &) const = 0;

  //! Query to see if data is in tabular or functional form (false).
  bool data_in_tabular_form() const { return false; }

  //! Query to get the reaction type.
  rtt_cdi::Reaction getReactionType() const { return reaction; }

  //! Query for model type.
  rtt_cdi::Model getModelType() const { return model; }

  //! Return the energy policy descriptor (mg).
  inline std_string getEnergyPolicyDescriptor() const {
    return std_string("mg");
  }

  // Get the data description of the opacity.
  virtual std_string getDataDescriptor() const = 0;

  // Get the name of the associated data file.
  inline std_string getDataFilename() const { return std_string(); }

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

  // Pack the Analytic_MultigroupOpacity into a character string.
  virtual sf_char pack() const = 0;

  /*!
   * \brief Returns the general opacity model type, defined in OpacityCommon.hh
   *
   * Since this is an analytic model, return 1 (rtt_cdi::ANALYTIC_TYPE)
   */
  rtt_cdi::OpacityModelType getOpacityModelType() const {
    return rtt_cdi::ANALYTIC_TYPE;
  }
};

} // end namespace rtt_cdi_analytic

#endif // __cdi_analytic_Analytic_MultigroupOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi_analytic/Analytic_MultigroupOpacity.hh
//---------------------------------------------------------------------------//
