//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/test/TestSolver.hh
 * \author Randy M. Roberts
 * \date   Wed Jan 26 12:58:31 2000
 * \brief  This class tests the rtt_LAMG::Solver class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMG_test_TestSolver_hh__
#define __LAMG_test_TestSolver_hh__

#include "Tester.hh"

// Forward References

namespace rtt_LAMG
{
class CompressedRowStorage;
class DynamicCompRowStorage;
} // end namespace rtt_LAMG

namespace rtt_LAMG_test
{
 
//===========================================================================//
/*!
 * \class TestSolver
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestSolver : public Tester
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestSolver(const std::string &filename, int &argc, char **argv,
	       std::ostream &os_in)
	: Tester(filename, argc, argv, os_in)
    {
	/* empty */
    }
    
    //DEFAULTED: TestSolver(const TestSolver &rhs);
    //DEFAULTED: ~TestSolver();

    // MANIPULATORS
    
    //DEFAULTED: TestSolver& operator=(const TestSolver &rhs);

    //! Main interface to testing class.
    
    void runTest();

    // ACCESSORS

    const std::string name() const { return "TestSolver"; }
	
  private:
    
    // IMPLEMENTATION

    rtt_LAMG::DynamicCompRowStorage makeDCRS(const int size);
    rtt_LAMG::CompressedRowStorage makeCRS(const int size);
};

} // end namespace rtt_LAMG_test

#endif                          // __LAMG_test_TestSolver_hh__

//---------------------------------------------------------------------------//
//                              end of LAMG/test/TestSolver.hh
//---------------------------------------------------------------------------//
