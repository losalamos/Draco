//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/tCDIGandolf.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:26:24 2000
 * \brief  Header file for tCDIGandolf.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_tCDIGandolf_hh__
#define __cdi_gandolf_tCDIGandolf_hh__

#include "UnitTestFrame/TestApp.hh"
#include "ds++/SP.hh"
#include "cdi/Opacity.hh"

#include <vector>

namespace rtt_cdi_gandolf_test
{
 
using std::vector;

//===========================================================================//
/*!
 * \class tCDIGandolf
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

class tCDIGandolf : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    tCDIGandolf( int argc, char *argv[], std::ostream &os_in );
    //Defaulted: tCDI(const tCDI &rhs);
    //Defaulted: ~tCDI();

    // MANIPULATORS
    
    //Defaulted: tCDI& operator=(const tCDI &rhs);

    // ACCESSORS
    std::string name() const { return "tCDIGandolf"; }
    std::string version() const;

  protected:
    
    std::string runTest();

  private:
    
    // IMPLEMENTATION
    
    /*!
     * \brief Returns true if the elements of the two vectors are
     *        identical to 10 decimal places.
     */
    bool match( const vector<double> computedValue, 
		const vector<double> refereneceValue );
    
    /*!
     * \brief Returns true if the two values are identical to 10
     *        decimal places. 
     */
    bool match( const double computedValue,
		const double referenceValue );

    // The CDIGandolf tests have been broken down in to smaller test
    // routines.  The prototypes for these routines are listed here.
    void testEnergyBoundaryAccessor( 
	const rtt_dsxx::SP<rtt_cdi::Opacity> spOpacity ); 
    void testDensityGridAccessor( 
	const rtt_dsxx::SP<rtt_cdi::Opacity> spOpacity ); 
    void testTemperatureGridAccessor( 
	const rtt_dsxx::SP<rtt_cdi::Opacity> spOpacity ); 
    template < class temperatureType, class densityType >
    bool testMGPlankOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi::Opacity> spOpacity,
	const temperatureType temperature,
	const densityType density,
	const vector<double> tabulatedValues,
	const std::string skey = "pmg" );
    template < class T1, class T2, class T3 >
    bool testGrayPlankOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi::Opacity> spOpacity,
	const T1 temperature, const T2 density,
	const T3 tabulatedValue, const std::string skey = "pgray" ); 
    template < class temperatureType, class densityType >
    bool testMGRosselandOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi::Opacity> spOpacity,
	const temperatureType temperature, const densityType density, 
	const vector<double> tabulatedValues,
	const std::string skey = "rtmg" ); 
    template < class T1, class T2, class T3 >
    bool testGrayRosselandOpacityAccessorPassed(
	const rtt_dsxx::SP<rtt_cdi::Opacity> spOpacity,
	const T1 temperature, const T2 density, 
	const T3 tabulatedValue, const std::string skey = "rgray" ); 
};

} // end namespace rtt_cdi_gandolf_test

#endif // __cdi_gandolf_tCDIGandolf_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/tCDI.hh
//---------------------------------------------------------------------------//
