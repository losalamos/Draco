//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/tCDI.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:26:24 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_tCDI_hh__
#define __cdi_tCDI_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_CDI_test
{
 
//===========================================================================//
/*!
 * \class tCDI
 *
 * \brief A class used to test the QuadCreator and Quadrature
 * classes.  It is inhierited from Randy's UnitTestFrame::TestApp.
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
};

    // ==================================================
    // ==================================================
    // Implement a dummy MatProps Class here???
    // ==================================================
    // ==================================================

} // end namespace rtt_cdi_test

#endif                          // __cdi_tCDI_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/tCDI.hh
//---------------------------------------------------------------------------//
