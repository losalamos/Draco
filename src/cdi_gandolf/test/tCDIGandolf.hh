//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/tCDIGandolf.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:26:24 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_tCDIGandolf_hh__
#define __cdi_gandolf_tCDIGandolf_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_cdi_gandolf_test
{
 
//===========================================================================//
/*!
 * \class tCDIGandolf
 *
 * \brief A class used to test the Gandolf Opacity reader under CDI
 * classes.  It is inhierited from Randy's UnitTestFrame::TestApp.
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
};

    // ==================================================
    // ==================================================
    // Implement a dummy MatProps Class here???
    // ==================================================
    // ==================================================

} // end namespace rtt_cdi_gandolf_test

#endif // __cdi_gandolf_tCDIGandolf_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/tCDI.hh
//---------------------------------------------------------------------------//
