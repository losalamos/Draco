//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/tGandolfWithCDI.hh
 * \author Kelly Thompson
 * \date   Sat Mar 10 16:02:27 2001
 * \brief  Header file for tGandolfWithCDI
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_tGandolfWithCDI_hh__
#define __cdi_tGandolfWithCDI_hh__

#include "UnitTestFrame/TestAppNoC4.hh"

namespace rtt_gandolf_with_cdi_test
{
    
    class tGandolfWithCDI : public rtt_UnitTestFrame::TestApp
    {
	
      public:
	
	// CREATORS
	
	tGandolfWithCDI( int argc, char *argv[], std::ostream &os_in );
	
	// ACCESSORS
	std::string name() const { return "tGandolfWithCDI"; }
	std::string version() const;
	
      protected:
	
	std::string runTest();
	
      private:
	
	// These routines are used by the test routine to compare computed 
	// results to tabulated values.  
	bool match( 
	    const double computedValue,
	    const double referenceValue ) const;
	bool match( 
	    const std::vector< double >& computedValue, 
	    const std::vector< double >& referenceValue ) const;
	bool match( 
	    const std::vector< std::vector< double > >& computedValue, 
	    const std::vector< std::vector< double > >& referenceValue ) const; 
    };
    
} // end namespace gandolf_with_cdi_test

#endif // __cdi_tGandolfWithCDI_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/tGandolfWithCDI.hh
//---------------------------------------------------------------------------//
