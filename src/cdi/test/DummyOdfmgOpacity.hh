//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOdfmgOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 17:12:51 2001
 * \brief  DummyOdfmgOpacity class header file (derived from
 *         ../OdfmgOpacity)
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.  */
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyOdfmgOpacity_hh__
#define __cdi_DummyOdfmgOpacity_hh__

#include "cdi/OdfmgOpacity.hh"

namespace rtt_cdi_test {

//===========================================================================//
/*!
 * \class DummyOdfmgOpacity
 *
 * \brief This is an opacity class that derives its interface from
 *        cdi/OdfmgOpacity and is used for testing purposes only.
 *
 * \sa This opacity class always contains the same data (set by the
 * default constructor).  The data table has the following properties:
 *
 *     Temperatures     = { 1.0, 2.0, 3.0 }
 *     Densities        = { 0.1, 0.2 }
 *     EnergyBoundaries = { 0.05, 0.5, 5.0, 50.0 }
 *     BandBoundaries   = { 0.00, 0.125, 0.25, 0.50, 1.00 }
 *
 *     Opacity = 2 * ( temperature + density/1000 )
 *                 / ( E_high + E_low )
 *
 * In addition to providing definitions for the member functions
 * outlined in OdfmgOpacity this class provides three additional
 * 1-D STL-like accessors for opacity data.
 *
 * \sa cdi/test/tDummyOpacity.cc
 * \sa cdi/test/tCDI.cc
 */
//===========================================================================//

class DLL_PUBLIC_cdi_test DummyOdfmgOpacity : public rtt_cdi::OdfmgOpacity {

  // DATA - all of these values are set in the constructor.

  // string descriptors
  std::string const dataFilename;           // "none"
  std::string const dataDescriptor;         // "DummyOdfmgOpacity"
  std::string const energyPolicyDescriptor; // "Odfmg"

  size_t const numTemperatures;    // = 3
  size_t const numDensities;       // = 2
  size_t const numGroupBoundaries; // = 4
  size_t const numBandBoundaries;  // = 5

  std::vector<double> groupBoundaries; // = { 0.05, 0.5, 5.0, 50.0 }
  std::vector<double> bandBoundaries;  // = { 0.00, 0.125, 0.25, 0.50, 1.00 }
  std::vector<double> temperatureGrid; // = { 1.0, 2.0, 3.0 }
  std::vector<double> densityGrid;     // = { 0.1, 0.2 }

  rtt_cdi::Reaction reaction_type;
  rtt_cdi::Model model_type;

public:
  // -------------------------- //
  // Constructors & Destrcutors //
  // -------------------------- //

  /*!
   * \brief Constructor for DummyOdfmgOpacity object.
   *
   * \sa The constructor assigns fixed values for all of the member data.
   *     Every instance of this object has the same member data.
   */
  DummyOdfmgOpacity(rtt_cdi::Reaction = rtt_cdi::TOTAL,
                    rtt_cdi::Model = rtt_cdi::ANALYTIC);

  /*!
   * \brief Constructor for DummyOdfmgOpacity object.
   *
   * This constructor allows the user to enter a different number of
   * frequency boundaries.
   */
  DummyOdfmgOpacity(rtt_cdi::Reaction reaction, rtt_cdi::Model model,
                    size_t num_groupboundaries, size_t num_bandboundaries);

  /*!
   * \brief Default DummyOdfmgOpacity destructor.
   *
   * This is required to correctly release memory when a DummyOdfmgOpacity
   * object is destroyed.
   */
  ~DummyOdfmgOpacity(void){/*empty*/};

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief Opacity accessor that returns a vector of opacities (one for
   *     each group) that corresponds to the provided temperature and
   *     density.
   *
   *     Opacity = 2 * ( temperature + density/1000 )
   *                 / ( E_high + E_low )
   *
   * \param targetTemperature The temperature value for which an opacity
   *     value is being requested (keV).
   * \param targetDensity The density value for which an opacity value is
   *     being requested (g/cm^3).
   * \return A vector of interpolated opacities (cm^2/g).  Each entry of
   *     this vector corresponds to one energy group.
   */
  std::vector<std::vector<double>> getOpacity(double targetTemperature,
                                              double targetDensity) const;

  /*!
   * \brief Opacity accessor that returns a vector of multigroup opacities
   *     corresponding to the provided vector of temperatures and a single
   *     density.  Each multigroup opacity is in itself a vector of
   *     numGroups opacities.
   *
   *     Opacity = 2 * ( temperature + density/1000 )
   *                 / ( E_high + E_low )
   *
   * \param targetTemperature A vector of temperature values for which
   *     corresponding opacity values are being requested (keV).
   * \param targetDensity The density value for which opacity values are
   *     being requested (g/cm^3).
   * \return A vector of interpolated multigroup opacities (cm^2/g).  Each
   *     entry of this vector corresponds to one of the provided
   *     temperatures and is itself a vector of numGroups opacities.
   */
  std::vector<std::vector<std::vector<double>>>
  getOpacity(const std::vector<double> &targetTemperature,
             double targetDensity) const;

  /*!
   * \brief Opacity accessor that returns a vector of multigroup opacities
   *     corresponding to the provided vector of densities and a single
   *     temperature.  Each multigroup opacity is in itself a vector of
   *     numGroups opacities.
   *
   *     Opacity = 2 * ( temperature + density/1000 )
   *                 / ( E_high + E_low )
   *
   * \param targetTemperature The temperature value for which corresponding
   *     opacity values are being requested (keV).
   *
   * \param targetDensity A vector of density values for which corresponding
   *     opacity values are being requested (g/cm^3).
   *
   * \return A vector of interpolated multigroup opacities (cm^2/g).  Each
   *     entry of this vector corresponds to one of the provided densities
   *     and is itself a vector of numGroups opacities.
   */
  std::vector<std::vector<std::vector<double>>>
  getOpacity(double targetTemperature,
             const std::vector<double> &targetDensity) const;

  /*!
   * \brief Opacity accessor that returns an STL container of opacities that
   *     correspond to a tuple of provided STL containers (temperatures and
   *     densities).  The length of the temperature and the the density
   *     container should be equal and the length of the opacity container
   *     should be numGroups x temperature.size().
   *
   *     This function is not required by OdfmgOpacity.
   *
   * \param tempIter The beginning position of a STL container that holds a
   *     list of temperatures (keV).
   * \param templast The end position of a STL container that holds a list
   *     of temperatures (keV).
   * \param densIter The beginning position of a STL container that holds a
   *     list of densities (g/cm^3).
   * \param densLast The end position of a STL container that holds a list
   *     of densities (g/cm^3).
   * \param opacityIter The beginning position of a STL container into
   *     which multigroup opacity values corresponding to the given tuple of
   *     (temperature, density) values and the number of energy groups will
   *     be stored (cm^2/g).
   * \return A list (of type OpacityIterator) of multigroup opacities are
   *     returned (cm^2/g).  These multigroup opacities correspond to the
   *     provided tuple of (temperature, density) values and the total
   *     number of energy groups.
   */
  template <class TemperatureIterator, class DensityIterator,
            class OpacityIterator>
  OpacityIterator getOpacity(TemperatureIterator tempIter,
                             TemperatureIterator templast,
                             DensityIterator densIter, DensityIterator densLast,
                             OpacityIterator opacityIter) const;

  /*!
   * \brief Opacity accessor that returns an STL container of opacities that
   *     correspond to a list of provided STL temperature values.  The
   *     length of the opacity container should be numGroups x
   *     temperature.size().
   *
   *     This function is not required by OdfmgOpacity.
   *
   * \param tempIter The beginning position of a STL container that holds a
   *     list of temperatures (keV).
   *
   * \param templast The end position of a STL container that holds a list
   *     of temperatures (keV).
   *
   * \param targetDensity A single density value for which opacity values
   *     are being requested (g/cm^3).
   *
   * \param opacityIter The beginning position of a STL container into
   *     which multigroup opacity values corresponding to the given list of
   *     temperature values and the number of energy groups will be stored
   *     (cm^2/g).
   *
   * \return A list (of type OpacityIterator) of multigroup opacities are
   *     returned (cm^2/g).  These multigroup opacities correspond to the
   *     provided list of temperature values, the fixed density and the
   *     total number of energy groups.
   */
  template <class TemperatureIterator, class OpacityIterator>
  OpacityIterator getOpacity(TemperatureIterator tempIter,
                             TemperatureIterator templast, double targetDensity,
                             OpacityIterator opacityIter) const;

  /*!
   * \brief Opacity accessor that returns an STL container of opacities that
   *     correspond to a list of provided STL density values and a fixed
   *     temperature.  The length of the opacity container should be
   *     numGroups x density.size().
   *
   *     This function is not required by OdfmgOpacity.
   *
   * \param targetTemperature A single temperature value for which opacity
   *     values are being requested (keV).
   * \param densIter The beginning position of a STL container that holds a
   *     list of densities (g/cm^3).
   * \param densLast The end position of a STL container that holds a list
   *     of densities (g/cm^3).
   * \param opacityIter The beginning position of a STL container into
   *     which multigroup opacity values corresponding to the given list of
   *     density values and the number of energy groups will be stored
   *     (cm^2/g).
   * \return A list (of type OpacityIterator) of multigroup opacities are
   *     returned (cm^2/g).  These multigroup opacities correspond to the
   *     provided list of density values, the fixed temperature and the
   *     total number of energy groups.
   */
  template <class DensityIterator, class OpacityIterator>
  OpacityIterator getOpacity(double targetTemperature, DensityIterator densIter,
                             DensityIterator densLast,
                             OpacityIterator opacityIter) const;

  /*!
   * \brief Data is in tables.
   */
  bool data_in_tabular_form() const { return true; }

  /*!
   * \brief Return the reaction type.
   */
  rtt_cdi::Reaction getReactionType() const { return reaction_type; }

  /*!
   * \brief Return the model type.
   */
  rtt_cdi::Model getModelType() const { return model_type; }

  /*!
   * \brief Returns a "plain English" description of the data.
   */
  std::string getDataDescriptor() const { return dataDescriptor; };

  /*!
   * \brief Returns a "plain English" description of the energy group
   *     structure (gray vs. multigroup).
   */
  std::string getEnergyPolicyDescriptor() const {
    return energyPolicyDescriptor;
  };

  /*!
   * \brief Returns the name of the associated data file.  Since there is no
   *     data file associated with this opacity class the string "none" is
   *     returned.
   */
  std::string getDataFilename() const { return dataFilename; };

  /*!
   * \brief Returns a vector of temperatures that define the cached opacity
   *     data table.
   */
  std::vector<double> getTemperatureGrid() const { return temperatureGrid; };

  /*!
   * \brief Returns a vector of densities that define the cached opacity
   *     data table.
   */
  std::vector<double> getDensityGrid() const { return densityGrid; };

  /*!
   * \brief Returns a vector of energy group boundaries that define the
   *     cached multigroup opacity data table.
   */
  std::vector<double> getGroupBoundaries() const { return groupBoundaries; };

  /*!
   * \brief Returns the size of the temperature grid.
   */
  size_t getNumTemperatures() const { return numTemperatures; };

  /*!
   * \brief Returns the size of the density grid.
   */
  size_t getNumDensities() const { return numDensities; };

  /*!
   * \brief Returns the number of energy group boundaries.
   */
  size_t getNumGroupBoundaries() const { return numGroupBoundaries; };

  /*!
   * \brief Returns the number of energy groups.
   */
  size_t getNumGroups() const { return numGroupBoundaries - 1; };

  /*!
   * \brief Returns a vector of points along the cumulative opacity
   *     distribution that mark the fraction of each band. First point is
   *     always zero, last point is always one.
   */
  std::vector<double> getBandBoundaries() const { return bandBoundaries; };

  /*!
   * \brief Returns the number of group boundaries found in the current
   *     multigroup data set.
   */
  size_t getNumBandBoundaries() const { return bandBoundaries.size(); };

  /*!
   * \brief Returns the number of band boundaries set in the current
   *     multigroup data set. The lower boundary is always zero, and the
   *     upper boundary is always one. Analagous to multigroup: numBands =
   *     numBandBoundaries - 1.
   */
  size_t getNumBands() const { return getNumBandBoundaries() - 1; };

  // Dummy pack function.
  std::vector<char> pack() const { return std::vector<char>(); }

  /*!
   * \brief Returns the general opacity model type, defined in
   *     OpacityCommon.hh
   */
  rtt_cdi::OpacityModelType getOpacityModelType() const {
    return rtt_cdi::DUMMY_TYPE;
  }

}; // end of class DummyOdfmgOpacity

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
// (enable us to use automatic instantiation)
//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns an STL container of opacities that
 *     correspond to a tuple of provided STL containers (temperatures and
 *     densities).  The length of the temperature and the the density
 *     container should be equal and the length of the opacity container
 *     should be numGroups x temperature.size().
 */
template <class TemperatureIterator, class DensityIterator,
          class OpacityIterator>
OpacityIterator DummyOdfmgOpacity::getOpacity(
    TemperatureIterator tempIter, TemperatureIterator tempLast,
    DensityIterator densIter, DensityIterator densLast,
    OpacityIterator opacityIter) const {
  size_t ng = numGroupBoundaries - 1;
  // loop over all temperatures and densities in the range
  // (tempFirst,tempLast) & (densIter,densLast).
  for (; densIter != densLast && tempIter != tempLast; ++tempIter, ++densIter)
    for (size_t ig = 0; ig < ng; ++ig, ++opacityIter)
      *opacityIter = 2.0 * (*tempIter + *densIter / 1000.0) /
                     (groupBoundaries[ig] + groupBoundaries[ig + 1]);
  return opacityIter;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns an STL container of opacities that
 *     correspond to a list of provided STL temperature values.  The length of
 *     the opacity container should be numGroups x temperature.size().
 */
template <class TemperatureIterator, class OpacityIterator>
OpacityIterator DummyOdfmgOpacity::getOpacity(
    TemperatureIterator tempIter, TemperatureIterator templast,
    double targetDensity, OpacityIterator opacityIter) const {
  size_t ng = numGroupBoundaries - 1;
  // loop over all temperatures in the range
  // (tempFirst,tempLast).
  for (; tempIter != templast; ++tempIter)
    for (size_t ig = 0; ig < ng; ++ig, ++opacityIter)
      *opacityIter = 2.0 * (*tempIter + targetDensity / 1000.0) /
                     (groupBoundaries[ig] + groupBoundaries[ig + 1]);
  return opacityIter;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns an STL container of opacities that
 *     correspond to a list of provided STL density values and a fixed
 *     temperature.  The length of the opacity container should be numGroups x
 *     density.size().
 */
template <class DensityIterator, class OpacityIterator>
OpacityIterator DummyOdfmgOpacity::getOpacity(
    double targetTemperature, DensityIterator densIter,
    DensityIterator densLast, OpacityIterator opacityIter) const {
  size_t ng = numGroupBoundaries - 1;
  // loop over all densities in the range
  // (densIter,densLast).
  for (; densIter != densLast; ++densIter)
    for (size_t ig = 0; ig < ng; ++ig, ++opacityIter)
      *opacityIter = 2.0 * (targetTemperature + *densIter / 1000.0) /
                     (groupBoundaries[ig] + groupBoundaries[ig + 1]);
  return opacityIter;
}

} // end namespace rtt_cdi_test

#endif // __cdi_DummyOdfmgOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/DummyOdfmgOpacity.hh
//---------------------------------------------------------------------------//
