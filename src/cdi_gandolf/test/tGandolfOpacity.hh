//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/tGandolfOpacity.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:26:24 2000
 * \brief  Header file for tCDIGandolf.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_tGandolfOpacity_hh__
#define __cdi_gandolf_tGandolfOpacity_hh__

#include "UnitTestFrame/TestApp.hh"
#include "ds++/SP.hh"

#include <vector>

namespace rtt_cdi_gandolf_test
{

//===========================================================================//
/*!
 * \class tGandolfOpacity
 *
 * \brief A class used to test the Gandolf Opacity reader under CDI
 * classes.  It is inherited from Randy's UnitTestFrame::TestApp.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class tGandolfOpacity : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    tGandolfOpacity( int argc, char *argv[], std::ostream &os_in );
    //Defaulted: tCDI(const tCDI &rhs);
    //Defaulted: ~tCDI();

    // MANIPULATORS
    
    //Defaulted: tCDI& operator=(const tCDI &rhs);

    // ACCESSORS
    std::string name() const { return "tGandolfOpacity"; }
    std::string version() const;

  protected:
    
    std::string runTest();

  private:
    
    // IMPLEMENTATION
    
    /*!
     * \brief Returns true if the two values are identical to 10
     *        decimal places. 
     */
    bool match( const double computedValue,
		const double referenceValue ) const;

    /*!
     * \brief Returns true if the elements of the two vectors are
     *        identical to 10 decimal places.
     */
    bool match( const std::vector< double >& computedValue, 
		const std::vector< double >& refereneceValue ) const;

    /*!
     * \brief Returns true if the two values are identical to 10
     *        decimal places. 
     */
    bool match( const std::vector< std::vector< double > >& computedValue,
 		const std::vector< std::vector< double > >& referenceValue ) const;

    
    // The CDIGandolf tests have been broken down in to smaller test
    // routines.  The prototypes for these routines are listed
    // here.
    
   /*!
     * \brief Test the getGroupBoundaries() accessor.
     *
     * \sa This test routine operates on the smart pointer to an
     *     opacity object provided as an argument.  It assumes that
     *     the OpacityObject has a Energy Group Structure with 13
     *     entries:  { 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1.0, 3.0, 7.0
     *                 10.0, 30.0, 70.0, 100.0 } keV
     */    
    template< class opacityClassType >
    void testEnergyBoundaryAccessor( const opacityClassType spOpacity ); 

    /*!
     * \brief Test the getDensityGrid() accessor.
     *
     * \sa This test routine operates on the smart pointer to an
     *     opacity object provided as an argument.  It assumes that
     *     the OpacityObject has a Density Grid with three entries 
     *     { 0.1, 0.5, 1.0 } g/cm^3.
     */
    template< class opacityClassType >
    void testDensityGridAccessor( const opacityClassType spOpacity ); 

    /*!
     * \brief Test the getTemperatureGrid() accessor.
     *
     * \sa This test routine operates on the smart pointer to an
     *     opacity object provided as an argument.  It assumes that
     *     the OpacityObject has a Temperature Grid with three entries 
     *     { 0.1, 1.0, 10.0 } keV.
     */
    template< class opacityClassType >
    void testTemperatureGridAccessor( const opacityClassType spOpacity ); 

    /*!
     * \brief Tests the getOpacity() function.
     *
     * \sa This templated function will use the getOpacity()
     *     accessor for the class specified by spOpacity.  The opacity
     *     values are interpolated using the information provided by
     *     "temperature" and "density."  The value(s) returned from
     *     the getOpacity() function is (are) assumed to be equal to
     *     the type of "tabulatedValues."  The contents of "tabulated
     *     Values" is then compared to the result of getOpacity().  If
     *     the two values match then the test returns "true" otherwise
     *     it returns "false."
     *
     * \param spOpacity A smart pointer to an Opacity object
     *     (i.e.: GandolfGrayOpacity or GandolfMultigroupOpacity) 
     *
     * \param temperature A list of temperatures (or a single value).
     *
     * \param density A list of densities (or a single value).
     *
     * \param tabulatedValues A list of opacity values (or a single
     *     value). that represent the exact values that should be
     *     returned by getOpacity.  The type of "tabulatedValues" must 
     *     match the type returned by getOpacity().
     */
    template < class temperatureType, class densityType, 
	class testValueType, class opacityClassType >
    bool opacityAccessorPassed(
 	const opacityClassType spOpacity,
 	const temperatureType temperature, 
	const densityType density, 
 	const testValueType tabulatedValues );
    
};

} // end namespace rtt_cdi_gandolf_test

#endif // __cdi_gandolf_tGandolfOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/tGandolfOpacity.hh
//---------------------------------------------------------------------------//
