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

#include <vector>
#include <string>

#include "../GrayOpacity.hh"

namespace rtt_dummyGrayOpacity
{

    /*!
     * \class DummyGrayOpacity
     *
     * \brief
     *
     * \sa This is a dummy opacity class that has the following
     *     properties: 
     *
     *     Temperatures = { 1.0, 2.0, 3.0 }
     *     Densities    = { 0.1, 0.2 }
     *     Opacity = temp + density/1000
     *
     * \example tDummyOpacity.cc
     * \example tCDI.cc
     */

    class DummyGrayOpacity : public rtt_cdiGrayOpacity::GrayOpacity
    {
	// DATA
	
 	const std::string dataFilename;
	const std::string dataDescriptor;
	const std::string energyPolicyDescriptor;
	
	const int numTemperatures;
	const int numDensities;
	
	std::vector< double > temperatureGrid;
	std::vector< double > densityGrid;

      public:

	// ------------ //
	// Constructors //
	// ------------ //
	
	/*!
	 * \brief Constructor for DummyOpacity object.
	 * 
	 * See DummyOpacity.hh for details.
	 *
	 * Note that everything in this file must be templated by the
	 * EnergyPolicy.  All Templated forms of DummyOpacity<EnergyPolicy>
	 * must be instantiated in DummyOpacity_pt.cc
	 */
	DummyGrayOpacity();
	 
	/*!
	 * \brief Default DummyOpacity() destructor.
	 *
	 * This is required to correctly release memory when a
	 * DummyOpacity<EnergyPolicy> is destroyed.
	 */
	~DummyGrayOpacity() {};

	// --------- //
	// Accessors //
	// --------- //
	
	/*!
	 * \brief Opacity accessor that returns a single opacity that 
	 *     corresponds to the provided temperature and density.
	 */
	double getOpacity( const double targetTemperature,
			   const double targetDensity ) const; 
	
	std::vector< double > getOpacity(
	    const std::vector< double >& targetTemperature,
	    const double targetDensity ) const;
	
	std::vector< double > getOpacity( 
	    const double targetTemperature,
	    const std::vector< double >& targetDensity ) const; 

	template< class OpacityIterator, class TemperatureIterator >
	OpacityIterator getOpacity( TemperatureIterator tempFirst,
				    TemperatureIterator tempLast,
				    const double targetDensity,
				    OpacityIterator opacityFirst ) const;

	template< class OpacityIterator, class DensityIterator >
	OpacityIterator getOpacity( const double targetTemperature,
				    DensityIterator densFirst,
				    DensityIterator densLast,
				    OpacityIterator opacityFirst ) const;

	// Given a (temperature,density) tuple return an opacity.
	template< class OpacityIterator, class TemperatureIterator,
  	          class DensityIterator >
	OpacityIterator getOpacity( TemperatureIterator tempFirst,
				    TemperatureIterator tempLast,
				    DensityIterator densFirst,
				    DensityIterator densLast,
				    OpacityIterator opacityFirst ) const;

	/*!
	 * \brief Returns a "plain English" description of the data.
	 */
	const std::string& getDataDescriptor() const { 
	    return dataDescriptor; };

	/*!
	 * \brief Returns a "plain English" description of the energy
	 *	  group structure (gray vs. multigroup).
	 */	
	const std::string& getEnergyPolicyDescriptor() const {
	    return energyPolicyDescriptor; };

	/*!
	 * \brief Returns the name of the associated IPCRESS file.
	 *
	 * The definition of this function is not included here to prevent 
	 *     the inclusion of the GandolfFile.hh definitions within this 
	 *     header file.
	 */
	const std::string& getDataFilename() const {
	    return dataFilename; };

	/*!
	 * \brief Returns a vector of temperatures that define the cached
	 *     opacity data table.
	 * 
	 * We do not return a const reference because this function
	 * must construct this information from more fundamental tables.
	 */
	std::vector<double> getTemperatureGrid() const {
	    return temperatureGrid; };
	
	/*!
	 * \brief Returns a vector of densities that define the cached
	 *     opacity data table.
	 * 
	 * We do not return a const reference because this function
	 * must construct this information from more fundamental tables.
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
