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

#include "UnitTestFrame/TestAppNoC4.hh"

namespace rtt_dummy_opacity_test
{

class tDummyOpacity : public rtt_UnitTestFrame::TestApp
{

  public:

    // CREATORS
    
    tDummyOpacity( int argc, char *argv[], std::ostream &os_in );

    // ACCESSORS
    std::string name() const { return "tDummyOpacity"; }
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

} // end namespace rtt_dummy_opacity_test

#endif // __cdi_dummy_tDummyOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/tDummyOpacity.hh
//---------------------------------------------------------------------------//
