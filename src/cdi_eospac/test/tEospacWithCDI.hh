//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/test/tEospacWithCDI.hh
 * \author Kelly Thompson
 * \date   Mon Apr 19 10:33:54 2001
 * \brief  Header file for tEospacWithCDI.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_tEospacWithCDI_hh__
#define __cdi_eospac_tEospacWithCDI_hh__

#include "UnitTestFrame/TestAppNoC4.hh"

namespace rtt_cdi_eospac_test
{
    
    //===========================================================================//
    /*!
     * \class tEospacWithCDI
     *
     * \brief A class used to test the Eospac object under the cdi_eospac
     *        package.  It inherits from Randy's
     *        UnitTestFrame::TestApp. This test program tests the use
     *        of cdi_eospac from within CDI.
     *
     */
    // revision history:
    // -----------------
    // 0) original
    // 
    //===========================================================================//
    
    class tEospac : public rtt_UnitTestFrame::TestApp
    {
	
	// NESTED CLASSES AND TYPEDEFS
	
	// DATA
	
      public:
	
	// CREATORS
	
	tEospac( int argc, char *argv[], std::ostream &os_in );
	
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
    };
    
} // end namespace rtt_cdi_eospac_test

#endif // __cdi_eospac_tEospacWithCDI_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/test/tEospacWithCDI.hh
//---------------------------------------------------------------------------//

