//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/test/TestConjGrad.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 10:19:48 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_test_TestConjGrad_hh__
#define __ConjGrad_test_TestConjGrad_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_ConjGrad_test
{
 
//===========================================================================//
/*!
 * \class TestConjGrad
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestConjGrad  : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestConjGrad(int argc, char *argv[], std::ostream &os_in);
    //Defaulted: TestConjGrad(const TestConjGrad &rhs);
    //Defaulted: ~TestConjGrad();

    // MANIPULATORS
    
    //Defaulted: TestConjGrad& operator=(const TestConjGrad &rhs);

    // ACCESSORS

    std::string name() const { return "TestConjGrad"; }

    std::string version() const;

  protected:
    
    std::string runTest();
    
  private:

    void runTestConjGrad(bool jacobiPrecon);
    void runTestConjGradMatVec();
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGrad_test

#endif                          // __ConjGrad_test_TestConjGrad_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/test/TestConjGrad.hh
//---------------------------------------------------------------------------//
