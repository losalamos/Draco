//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_dummy/tDummyEoS.hh
 * \author Kelly Thompson
 * \date   Mon Apr 16 13:23:08 2001
 * \brief  Header file for tDummyEoS
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_tDummyEoS_hh__
#define __cdi_tDummyEoS_hh__

#include "UnitTestFrame/TestAppNoC4.hh"

namespace rtt_dummy_eos_test
{

    class tDummyEoS : public rtt_UnitTestFrame::TestApp
    {
	
      public:
	
	// CREATORS
	
	tDummyEoS( int argc, char *argv[], std::ostream &os_in );
	
	// ACCESSORS
	std::string name() const { return "tDummyEoS"; }
	std::string version() const;
	
      protected:
	
	std::string runTest();
	
      private:
	
	// These routines are used by the test routine to compare computed 
	// results to tabulated values.  
	bool match( 
	    double computedValue,
	    double referenceValue ) const;
	bool match( 
	    const std::vector< double >& computedValue, 
	    const std::vector< double >& referenceValue ) const;
    };
    
} // end namespace rtt_dummy_eos_test

#endif // __cdi_dummy_tDummyEoS_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/tDummyEoS.hh
//---------------------------------------------------------------------------//
