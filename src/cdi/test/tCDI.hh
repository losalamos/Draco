//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/tCDI.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:26:24 2000
 * \brief  Header file for the CDI class unit test
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_tCDI_hh__
#define __cdi_tCDI_hh__

#include "UnitTestFrame/TestAppNoC4.hh"

#include <vector>
#include <string>

namespace rtt_CDI_test
{

//===========================================================================//
/*!
 * \class tCDI
 *
 * \brief A class used to test the QuadCreator and Quadrature
 *        classes.  It is inhierited from Randy's
 *        UnitTestFrame::TestApp. 
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class tCDI : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    tCDI( int argc, char *argv[], std::ostream &os_in );

    // ACCESSORS
    std::string name() const { return "tCDI"; }
    std::string version() const;

  protected:
    
    std::string runTest();

  private:
    
    // IMPLEMENTATION

    // These routines are used by the test routine to compare computed 
    // results to tabulated values.  
    bool match( 
	const double computedValue,
	const double referenceValue ) const;
    bool match( 
	const std::vector< double >& computedValue, 
	const std::vector< double >& referenceValue ) const;
};

} // end namespace rtt_cdi_test

#endif // __cdi_tCDI_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/tCDI.hh
//---------------------------------------------------------------------------//
