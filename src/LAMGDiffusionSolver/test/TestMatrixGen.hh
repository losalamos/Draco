//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/TestMatrixGen.hh
 * \author Randy M. Roberts
 * \date   Fri Aug 20 09:11:42 1999
 * \brief  Header file for the TestMatrixGen class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_test_TestMatrixGen_hh__
#define __LAMGDiffusionSolver_test_TestMatrixGen_hh__

#include "Tester.hh"
#include "../MatrixP1DiffTraits.hh"
#include "../MatrixP1Diff.hh"

#include <iosfwd>
#include <string>

namespace rtt_LAMGDiffusionSolver_test
{
 
//===========================================================================//
/*!
 * \class TestMatrixGen
 *
 * \brief  This class tests whether a given MT is a model of the
 *            "Solon" MT concept.
 *
 * The TestMatrixGen class is templated on an MTFactory.
 * The MTFactory is used to create coupled instances of meshes and field
 * constructors, along with typedefs that determine whether the mesh is
 * structured or unstructured.
 *
 * Required MTFactory services include:
 *
 * MTFactory::MT -- The mesh type
 *
 * MTFactory::FieldConstructor -- The field constructor type.
 *
 * MTFactory::Product -- A class that is returned by the MTFactroy::create()
 * method.  This class is responsible for returning references to the mesh
 * and field constructor via the mesh() and fieldConstructor() methods,
 * respectively.
 *
 * MTFactory::Structured -- A tag that, if used, implies that the mesh
 * is a structured mesh.
 *
 * MTFactory::UnStructured -- A tag that, if used, implies that the mesh
 * is an unstructured mesh.
 *
 * MTFactory::Structuring -- A typedef to either, Structured, or UnStructured,
 * used to determine whether the mesh is actually one or the other.
 *
 * MTFactory::Product MTFactory::create() -- The method that returns a new
 * pair of meshes and field constructors.
 */
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MTFactory>
class TestMatrixGen : public Tester<MTFactory>
{

    // NESTED CLASSES AND TYPEDEFS

    typedef rtt_LAMGDiffusionSolver::MatrixP1Diff Matrix;
    typedef rtt_traits::MatrixFactoryTraits<Matrix> MatFacTraits;
    typedef typename MatFacTraits::PreComputedState PCMS;

    // DATA

  public:

    // CREATORS

    //! Constructor
    
    TestMatrixGen(const std::string &filename, int &argc, char **argv,
		   std::ostream &os_in)
	: Tester<MTFactory>(filename, argc, argv, os_in)
    {
	/* empty */
    }

    //! Destructor
 
    ~TestMatrixGen() { /* empty */ }

    // MANIPULATORS

    //! Main interface to testing class.
    
    void runTest();

    // ACCESSORS

    const std::string name() const { return "TestMatrixGen"; }
	
  private:
    
    // DISSALLOWED CREATORS

    TestMatrixGen(const TestMatrixGen &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMatrixGen& operator=(const TestMatrixGen &rhs);

    // IMPLEMENTATION

    void verifyPreComputedState(const PCMS &pcms);

    void verifyMultiplication(const std::vector<double> &b,
			      const typename TestMatrixGen::MT &mesh,
			      const typename TestMatrixGen::MT::ccsf &Adiag,
			      const typename TestMatrixGen::MT::fcdsf &Aoff,
			      const std::vector<double> &x);

    void setDiagonal(typename TestMatrixGen::MT::ccsf &val)
    {
	val = 6;
    }
    void setOffDiagonal(typename TestMatrixGen::MT::fcdsf &val)
    {
	val = -1;
    }
};

} // end namespace rtt_LAMGDiffusionSolver_test

#endif // __LAMGDiffusionSolver_test_TestMatrixGen_hh__

//---------------------------------------------------------------------------//
//                 end of LAMGDiffusionSolver/test/TestMatrixGen.hh
//---------------------------------------------------------------------------//
