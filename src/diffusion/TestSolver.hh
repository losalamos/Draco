//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diffusion/TestSolver.hh
 * \author Randy M. Roberts
 * \date   Mon Feb  7 10:26:24 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __diffusion_TestSolver_hh__
#define __diffusion_TestSolver_hh__

#include "Tester.hh"
#include "Diffusion_DB.hh"
#include "traits/MatrixFactoryTraits.hh"

namespace rtt_diffusion
{

// Forward Reference

template<class MT> class P1Matrix;

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

template<class MT, class Solver>
class TestSolver : public Tester
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::ccsf ccsf;
    typedef typename MT::fcdsf fcdsf;

    typedef typename Solver::Matrix Matrix;
    typedef rtt_traits::MatrixFactoryTraits<Matrix> MatFacTraits;
     

    typedef typename MT::FieldConstructor FieldConstructor;
    
    // DATA
    
    FieldConstructor fCtor;
    Solver &solver;
    Diffusion_DB diff_db;
    double tolerance;
    std::string version_m;
     
  public:

    // CREATORS
    
    TestSolver(const std::string &solverName, int &argc,
	       char *argv[], const std::string &version_in,
	       const FieldConstructor &fCtor_in, Solver &solver_in,
	       const Diffusion_DB &diff_db_in,
	       double tolerance_in, std::ostream &os_in)
	: Tester(std::string("TestSolver<") + solverName + ">", os_in,
		 argc, argv),
	  fCtor(fCtor_in), solver(solver_in), diff_db(diff_db_in),
	  tolerance(tolerance_in), version_m(version_in)
    {
	// empty
    }
    
    //DEFAULTED: TestSolver(const TestSolver &rhs);
    //DEFAULTED: ~TestSolver();

    // MANIPULATORS
    
    //DEFAULTED: TestSolver& operator=(const TestSolver &rhs);

    // ACCESSORS

  protected:
    
    // PROTECTED IMPLEMENTATION

    const std::string &version() const { return version_m; }
    
    void runTest();

  private:
    
    // IMPLEMENTATION

    void runTest(bool jacobiScale);
    P1Matrix<MT> createP1Matrix() const;
};

} // end namespace rtt_diffusion

#endif                          // __diffusion_TestSolver_hh__

//---------------------------------------------------------------------------//
//                              end of diffusion/TestSolver.hh
//---------------------------------------------------------------------------//
