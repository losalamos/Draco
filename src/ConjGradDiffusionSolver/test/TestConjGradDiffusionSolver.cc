//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGradDiffusionSolver/test/TestConjGradDiffusionSolver.cc
 * \author Randy M. Roberts
 * \date   Mon Apr 24 10:19:48 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestConjGradDiffusionSolver.hh"

#include "../SolverP1Diff.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"

#include "diffusion/TestSolver.t.hh"
#include "mesh/Mesh_XYZ.hh"
#include "nml/Group.hh"
#include "nml/Item.hh"
#include "nml/Items.hh"

#include <sstream>

namespace
{

class Test_DB {
  public:

#include ".nml_test.hh"

  public:
    void setup_namelist( NML_Group& g )
    {
#include ".nml_test.cc"
    }
};

template<class FT, class MAT>
inline void multiply(FT &b, const MAT &A, const FT &x)
{
    A.multiply(b, x);
}

} // end unnamed namespace

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_ConjGradDiffusionSolver_test::TestConjGradDiffusionSolver;
    
    return SP<TestApp>(new TestConjGradDiffusionSolver(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_ConjGradDiffusionSolver_test
{

using std::string;

TestConjGradDiffusionSolver::TestConjGradDiffusionSolver(int argc,
							 char *argv[],
							 std::ostream& os_in)
    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestConjGradDiffusionSolver" << std::endl;
}

string TestConjGradDiffusionSolver::version() const
{
    return rtt_ConjGradDiffusionSolver::release();
}

string TestConjGradDiffusionSolver::runTest()
{
    using rtt_diffusion::TestSolver;
    using rtt_ConjGradDiffusionSolver::SolverP1Diff;

    typedef TestSolver<Mesh_XYZ,SolverP1Diff<Mesh_XYZ> > Tester;

    NML_Group g("testConjGradDiffusionSolver");

    rtt_diffusion::Diffusion_DB diff_db;
    diff_db.setup_namelist(g);

    Test_DB tdb;
    tdb.setup_namelist(g);
    
    rtt_mesh::Mesh_DB mdb;
    mdb.setup_namelist(g);

    g.readgroup("testConjGradDiffusionSolver.in");
    g.writegroup("testConjGradDiffusionSolver.out");

    mdb.resize();
    
    rtt_dsxx::SP<Mesh_XYZ> spmesh(new Mesh_XYZ(mdb));

    const Mesh_XYZ::FieldConstructor &fCtor = spmesh;

    bool usePreconditioner = true;
    bool verbose = true;
    SolverP1Diff<Mesh_XYZ> solver(spmesh, fCtor, tdb.maxIterations,
				  tdb.epsilon, usePreconditioner, verbose);

    int argc = 0;
    char **argv = 0;
    
    Tester tester("ConjGradDiffusionSolver::SolverP1Diff", argc, argv,
		  rtt_ConjGradDiffusionSolver::release(),
		  fCtor, solver, diff_db, tdb.tolerance, os());
    tester.run();

    if (tester.passed())
	pass() << "Test ran to completion.";
    else
	fail() << "Test failed to run to completion.";

    if (passed())
    {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

} // end namespace rtt_ConjGradDiffusionSolver_test


//---------------------------------------------------------------------------//
//                              end of TestConjGradDiffusionSolver.cc
//---------------------------------------------------------------------------//
