//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyMultigroupOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 17:12:51 2001
 * \brief  DummyMultigroupOpacity class header file (derived from 
 *         cdi/MultigroupOpacity)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyMultigroupOpacity_hh__
#define __cdi_DummyMultigroupOpacity_hh__

#include <vector>
#include <string>

#include "../MultigroupOpacity.hh"

namespace rtt_dummyMultigroupOpacity
{

    /*!
     * \class DummyMultigroupOpacity
     *
     * \brief
     *
     * \sa This is a dummy opacity class that has the following
     *     properties: 
     *
     *     Temperatures = { 1.0, 2.0, 3.0 }
     *     Densities    = { 0.1, 0.2 }
     *     GroupBoundaries = { 0.05, 0.5, 5.0, 50.0 }
     *
     *     Opacity = 2.0 * ( Temperature + Density/1000.0 ) 
     *                   / ( E_high + E_low )
     *
     * \example tDummyOpacity.cc
     * \example tCDI.cc
     */

class DummyMultigroupOpacity : public rtt_cdiMultigroupOpacity::MultigroupOpacity
{

    // DATA
    
    const std::string dataFilename;
    const std::string dataDescriptor;
    const std::string energyPolicyDescriptor;

    const int numTemperatures;
    const int numDensities;
    const int numGroupBoundaries;

    std::vector< double> groupBoundaries;
    std::vector< double > temperatureGrid;
    std::vector< double > densityGrid;

  public:

    // ------------ //
    // Constructors //
    // ------------ //

    DummyMultigroupOpacity( );
    ~DummyMultigroupOpacity( ) { };

    // --------- //
    // Accessors //
    // --------- //

    const std::string& getDataDescriptor() const {
	return dataDescriptor; };

    template < class TemperatureIterator, class DensityIterator, 
	       class OpacityIterator >
    OpacityIterator getOpacity( TemperatureIterator tempIter,
				TemperatureIterator templast,
				DensityIterator densIter,
				DensityIterator densLast,
				OpacityIterator opacityIter ) const;
    
    template < class TemperatureIterator, class OpacityIterator >
    OpacityIterator getOpacity( TemperatureIterator tempIter,
				TemperatureIterator templast,
				const double targetDensity,
				OpacityIterator opacityIter ) const;

    template < class DensityIterator, class OpacityIterator >
    OpacityIterator getOpacity( const double targetTemperature,
				DensityIterator densIter,
				DensityIterator densLast,
				OpacityIterator opacityIter ) const;

    std::vector< double > getOpacity( const double targetTemperature,
				      const double targetDensity ) const; 

    std::vector< std::vector <double > > getOpacity( 
	const std::vector< double >& targetTemperature,
	const double targetDensity ) const; 

    std::vector< std::vector< double > > getOpacity( 
	const double targetTemperature,
	const std::vector< double >& targetDensity ) const; 

    const std::string& getEnergyPolicyDescriptor() const { 
	return energyPolicyDescriptor; };

    const std::string& getDataFilename() const {
	return dataFilename; };

    std::vector<double> getTemperatureGrid() const {
	return temperatureGrid; };

    std::vector<double> getDensityGrid() const {
	return densityGrid; };

    const std::vector<double>& getGroupBoundaries() const {
	return groupBoundaries; };
    
    int getNumTemperatures() const {
	return numTemperatures; };

    int getNumDensities() const {
	return numDensities; };

    int getNumGroupBoundaries() const {
	return numGroupBoundaries; };

}; // end of class GandolfOpacity

} // end namespace rtt_dummyMultigroupOpacity

#endif // __cdi_DummyMultigroupOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/DummyMultigroupOpacity.hh
//---------------------------------------------------------------------------//
