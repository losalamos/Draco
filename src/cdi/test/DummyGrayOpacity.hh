//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyGrayOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:29:17 2001
 * \brief  DummyGrayOpacity class header file (derived from ../GrayOpacity)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyGrayOpacity_hh__
#define __cdi_DummyGrayOpacity_hh__

#include "../GrayOpacity.hh"
#include "../OpacityCommon.hh"

namespace rtt_dummyGrayOpacity
{

//========================================================================
/*!
 * \class DummyGrayOpacity
 *
 * \breif This is an opacity class that derives its interface from
 * cdi/GrayOpacity and is used for testing purposes only.
 *
 * \sa This opacity class always contains the same data (set by the
 * default constructor).  The data table has the following properties:
 *
 *     Temperatures = { 1.0, 2.0, 3.0 }
 *     Densities    = { 0.1, 0.2 }
 *
 *     Opacity = temperature + density/1000
 *
 * In addition to providing definitions for the member functions
 * outlined in GrayOpacity this class provides three additional 1-D
 * STL-like accessors for opacity data.
 */
    
/*!
 * \example cdi/test/tDummyOpacity.cc
 * \example cdi/test/tCDI.cc
 */
//========================================================================

class DummyGrayOpacity : public rtt_cdi::GrayOpacity
{
    // DATA - all of these values are set in the constructor.
	
    // string descriptors
    const std::string dataFilename;            // "none"
    const std::string dataDescriptor;          // "DummyGrayOpacity"
    const std::string energyPolicyDescriptor;  // "Gray"
	
    // data grid size
    const int numTemperatures;  // = 3
    const int numDensities;     // = 2
	
    // the data grid
    std::vector< double > temperatureGrid;  // = { 1.0, 2.0, 3.0 }
    std::vector< double > densityGrid;      // = { 0.1, 0.2 }
	
  public:
	
    // -------------------------- //
    // Constructors & Destructors //
    // -------------------------- //
	
    /*!
     * \brief Constructor for DummyGrayOpacity object.
     * 
     * The constructor assigns fixed values for all of the member
     * data.  Every instance of this object has the same member
     * data. 
     */
    DummyGrayOpacity();
	
    /*!
     * \brief Default DummyGrayOpacity() destructor.
     *
     * This is required to correctly release memory when a
     * DummyGrayOpacity object is destroyed.
     */
    ~DummyGrayOpacity() {};
	
    // --------- //
    // Accessors //
    // --------- //
	
    /*!
     * \brief Opacity accessor that returns a single opacity that 
     *     corresponds to the provided temperature and density.  
     *
     *     Opacity = temperature + density/1000
     *
     * \param targetTemperature The temperature value for which an
     *     opacity value is being requested (keV).
     * \param targetDensity The density value for which an opacity 
     *     value is being requested (g/cm^3).
     * \return A single interpolated opacity (cm^2/g).
     */
    double getOpacity( double targetTemperature,
		       double targetDensity ) const; 
	
    /*!
     * \brief Opacity accessor that returns a vector of opacities that
     *     correspond to the provided vector of temperatures and a
     *     single density value. 
     *
     *     Opacity[i] = temperature[i] + density/1000
     *
     * \param targetTemperature A vector of temperature values for
     *     which opacity values are being requested (keV).
     * \param targetDensity The density value for which an opacity 
     *     value is being requested (g/cm^3).
     * \return A vector of opacities (cm^2/g).
     */
    std::vector< double > getOpacity(
	const std::vector< double >& targetTemperature,
	double targetDensity ) const;
	
    /*!
     * \brief Opacity accessor that returns a vector of opacities
     *     that correspond to the provided vector of densities and a
     *     single temperature value. 
     *
     *     Opacity[i] = temperature[i] + density/1000
     *
     * \param targetTemperature The temperature value for which an 
     *     opacity value is being requested (keV).
     * \param targetDensity A vector of density values for which
     *     opacity values are being requested (g/cm^3).
     * \return A vector of opacities (cm^2/g).
     */
    std::vector< double > getOpacity( 
	double targetTemperature,
	const std::vector< double >& targetDensity ) const; 
	
    /*! 
     * \brief Opacity accessor that returns an STL container of
     *     opacities that correspond to the provided STL container of
     *     temperatures and a single density.  The length of the
     *     opacity container and the temperature container should be
     *     equal. 
     *
     *     This function is not required by GrayOpacity.
     *
     * \param tempFirst The beginning position of a STL container
     *     that holds a list of temperatures (keV).
     * \param tempLast The end position of a STL container that
     *     holds a list of temperatures (keV).
     * \param targetDensity The density value for which an opacity 
     *     value is being requested (g/cm^3).
     * \param opacityFirst The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given temperature values will be stored (cm^2/g).
     * \return A list (of type OpacityIterator) of opacities are
     *     returned (cm^2/g).  These opacities correspond to the
     *     provided list of temperatures and the fixed density.
     */
    template< class OpacityIterator, class TemperatureIterator >
    OpacityIterator getOpacity( TemperatureIterator tempFirst,
				TemperatureIterator tempLast,
				const double targetDensity,
				OpacityIterator opacityFirst ) const;
    
    /*! 
     * \brief Opacity accessor that returns an STL container of
     *     opacities that correspond to the provided STL container of
     *     densities and a single temperature.  The length of the
     *     opacity container and the density container should be
     *     equal. 
     *
     *     This function is not required by GrayOpacity.
     *
     * \param targetTemperature The temperature value for which an
     *     opacity value is being requested (keV).
     * \param densFirst The beginning position of a STL container
     *     that holds a list of densities (g/cm^3).
     * \param densLast The end position of a STL container that
     *     holds a list of densities (g/cm^3).
     * \param opacityFirst The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given density values will be stored (cm^2/g).
     * \return A list (of type OpacityIterator) of opacities are
     *     returned (cm^2/g).  These opacities correspond to the
     *     provided list of densities and the fixed temperature.
     */
    template< class OpacityIterator, class DensityIterator >
    OpacityIterator getOpacity( double targetTemperature,
				DensityIterator densFirst,
				DensityIterator densLast,
				OpacityIterator opacityFirst ) const;
    
    /*! 
     * \brief Opacity accessor that returns an STL container of
     *     opacities that correspond to a tuple of provided STL
     *     containers (temperatures and densities).  The length of the
     *     opacity container, the temperature and the the density
     *     container should be equal. 
     *
     *     This function is not required by GrayOpacity.
     *
     * \param tempFirst The beginning position of a STL container
     *     that holds a list of temperatures (keV).
     * \param tempLast The end position of a STL container that
     *     holds a list of temperatures (keV).
     * \param densFirst The beginning position of a STL container
     *     that holds a list of densities (g/cm^3).
     * \param densLast The end position of a STL container that
     *     holds a list of densities (g/cm^3).
     * \param opacityFirst The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given tuple of (temperature, density) values will be stored
     *     (cm^2/g). 
     * \return A list (of type OpacityIterator) of opacities are
     *     returned (cm^2/g).  These opacities correspond to the
     *     provided tuple of (temperature, density) values.
     */
    template< class OpacityIterator, class TemperatureIterator,
	      class DensityIterator >
    OpacityIterator getOpacity( TemperatureIterator tempFirst,
				TemperatureIterator tempLast,
				DensityIterator densFirst,
				DensityIterator densLast,
				OpacityIterator opacityFirst ) const;

    /*!
     * \brief Data is analytic, not in tables.
     */
    bool data_in_tabular_form() const { return false; }

    /*!
     * \brief Return the reaction type.
     */
    rtt_cdi::Reaction getReactionType() const { return rtt_cdi::TOTAL; }

    /*!
     * \brief Returns a "plain English" description of the data.
     */
    std::string getDataDescriptor() const { 
	return dataDescriptor; };

    /*!
     * \brief Returns a "plain English" description of the energy
     *	  group structure (gray vs. multigroup).
     */	
    std::string getEnergyPolicyDescriptor() const {
	return energyPolicyDescriptor; };
	
    /*!
     * \brief Returns the name of the associated data file.  Since
     *     there is no data file associated with this opacity class
     *     the string "none" is returned.
     */
    std::string getDataFilename() const {
	return dataFilename; };
	
    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     */
    std::vector<double> getTemperatureGrid() const {
	return temperatureGrid; };
	
    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     */
    std::vector<double> getDensityGrid() const {
	return densityGrid; };
	
    /*!
     * \brief Returns the size of the temperature grid.
     */
    int getNumTemperatures() const { return numTemperatures; };
	
    /*! 
     * \brief Returns the size of the density grid.
     */
    int getNumDensities() const { return numDensities; };
	
};
    
} // end namespace rtt_dummyGrayOpacity

#endif // __cdi_DummyGrayOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/DummyGrayOpacity.hh
//---------------------------------------------------------------------------//
