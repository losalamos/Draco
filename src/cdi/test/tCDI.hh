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

#include "UnitTestFrame/TestApp.hh"

#include <vector>
#include <string>

namespace rtt_CDI_test
{

using std::vector;
using std::string;
 
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
    //Defaulted: tCDI(const tCDI &rhs);
    //Defaulted: ~tCDI();

    // MANIPULATORS
    
    //Defaulted: tCDI& operator=(const tCDI &rhs);

    // ACCESSORS
    std::string name() const { return "tCDI"; }
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

};

} // end namespace rtt_cdi_test

#endif                          // __cdi_tCDI_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/tCDI.hh
//---------------------------------------------------------------------------//
