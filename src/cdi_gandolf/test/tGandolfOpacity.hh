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
#include "cdi/Opacity.hh"
#include "../GandolfOpacity.hh"

#include <vector>

namespace rtt_cdi_gandolf_test
{
 
//     // forward declaration
//     class rtt_cdi_gandolf::GrayOpacity;
//     doesn't work!!!!!!!!!!!!
 

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
	const rtt_dsxx::SP<rtt_cdi_gandolf::MultigroupOpacity> spOpacity ); 
    void testDensityGridAccessor( 
	const rtt_dsxx::SP<rtt_cdi_gandolf::MultigroupOpacity> spOpacity ); 
    void testTemperatureGridAccessor( 
	const rtt_dsxx::SP<rtt_cdi_gandolf::MultigroupOpacity> spOpacity ); 
    template < class temperatureType, class densityType, class TestValueType >
    bool testMGPlankOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi_gandolf::MultigroupOpacity> spOpacity,
	const temperatureType temperature,
	const densityType density,
	const TestValueType tabulatedValues );
    template < class T1, class T2, class T3 >
    bool testGrayPlankOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi_gandolf::GrayOpacity> spOpacity,
	const T1 temperature, const T2 density,
	const T3 tabulatedValue ); 
    template < class temperatureType, class densityType, class TestValueType >
    bool testMGRosselandOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi_gandolf::MultigroupOpacity> spOpacity,
	const temperatureType temperature, const densityType density, 
	const TestValueType tabulatedValues );
    template < class T1, class T2, class T3 >
    bool testGrayRosselandOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi_gandolf::GrayOpacity> spOpacity,
	const T1 temperature, const T2 density, 
	const T3 tabulatedValue ); 

};

} // end namespace rtt_cdi_gandolf_test

#endif // __cdi_gandolf_tGandolfOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/tGandolfOpacity.hh
//---------------------------------------------------------------------------//
