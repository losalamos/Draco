//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_dummy/tDummyOpacity.hh
 * \author Kelly Thompson
 * \date   Fri Jan 5 16:00:24 2001
 * \brief  Header file for tDummyOpacity
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_tDummyOpacity_hh__
#define __cdi_tDummyOpacity_hh__

#include "UnitTestFrame/TestApp.hh"
#include "ds++/SP.hh"
#include "DummyOpacity.hh"

#include <vector>

namespace rtt_cdi_dummy_opacity_test
{
 
//===========================================================================//
/*!
 * \class tDummyOpacity
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

class tDummyOpacity : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    tDummyOpacity( int argc, char *argv[], std::ostream &os_in );
    //Defaulted: tCDI(const tCDI &rhs);
    //Defaulted: ~tCDI();

    // MANIPULATORS
    
    //Defaulted: tCDI& operator=(const tCDI &rhs);

    // ACCESSORS
    std::string name() const { return "tDummyOpacity"; }
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
		const double referenceValue ) const ;

    /*!
     * \brief Returns true if the elements of the two vectors are
     *        identical to 10 decimal places.
     */
    bool match( const std::vector<double>& computedValue, 
		const std::vector<double>& refereneceValue ) const;

    /*!
     * \brief Returns true if the two values are identical to 10
     *        decimal places. 
     */
    bool match( const std::vector< std::vector< double > >& computedValue,
		const std::vector< std::vector< double > >& referenceValue ) const;

    
    // The CDIGandolf tests have been broken down in to smaller test
    // routines.  The prototypes for these routines are listed here.
    void testEnergyBoundaryAccessor( 
	const rtt_dsxx::SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity ); 
    void testDensityGridAccessor( 
	const rtt_dsxx::SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity ); 
    void testTemperatureGridAccessor( 
	const rtt_dsxx::SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity ); 
    template < class temperatureType, class densityType, class TestValueType >
    bool testMGPlankOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity,
	const temperatureType temperature,
	const densityType density,
	const TestValueType tabulatedValues );
    template < class T1, class T2, class T3 >
    bool testGrayPlankOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_dummy_opacity::GrayOpacity> spOpacity,
	const T1 temperature, const T2 density,
	const T3 tabulatedValue ); 
    template < class temperatureType, class densityType, class TestValueType >
    bool testMGRosselandOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity,
	const temperatureType temperature, const densityType density, 
	const TestValueType tabulatedValues );
    template < class T1, class T2, class T3 >
    bool testGrayRosselandOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_dummy_opacity::GrayOpacity> spOpacity,
	const T1 temperature, const T2 density, 
	const T3 tabulatedValue ); 

};

} // end namespace rtt_dummy_opacity_test

#endif // __cdi_dummy_tDummyOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/tDummyOpacity.hh
//---------------------------------------------------------------------------//
