//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGradDiffusionSolver/test/TestConjGradDiffusionSolver.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 10:19:48 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGradDiffusionSolver_test_TestConjGradDiffusionSolver_hh__
#define __ConjGradDiffusionSolver_test_TestConjGradDiffusionSolver_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_ConjGradDiffusionSolver_test
{
 
//===========================================================================//
/*!
 * \class TestConjGradDiffusionSolver
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestConjGradDiffusionSolver  : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestConjGradDiffusionSolver(int argc, char *argv[], std::ostream &os_in);
    //Defaulted: TestConjGradDiffusionSolver(
    //                    const TestConjGradDiffusionSolver &rhs);
    //Defaulted: ~TestConjGradDiffusionSolver();

    // MANIPULATORS
    
    //Defaulted: TestConjGradDiffusionSolver&
    //                      operator=(const TestConjGradDiffusionSolver &rhs);

    // ACCESSORS

    std::string name() const { return "TestConjGradDiffusionSolver"; }

    std::string version() const;

  protected:
    
    std::string runTest();
    
  private:

    // IMPLEMENTATION
};

} // end namespace rtt_ConjGradDiffusionSolver_test

#endif // __ConjGradDiffusionSolver_test_TestConjGradDiffusionSolver_hh__

//---------------------------------------------------------------------------//
// end of ConjGradDiffusionSolver/test/TestConjGradDiffusionSolver.hh
//---------------------------------------------------------------------------//
