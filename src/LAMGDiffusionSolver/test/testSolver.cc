//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/testSolver.cc
 * \author Randy M. Roberts
 * \date   Mon Feb  7 10:49:55 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../SolverP1Diff.hh"
#include "../MatrixP1Diff.hh"
#include "../MatrixP1DiffTraits.hh"
#include "../Release.hh"
#include "MatrixP1DiffUtils.hh"

#include "diffusion/TestSolver.t.hh"
#include "mesh/Mesh_XYZ.hh"
#include "mesh/Mesh_XYZTraits.hh"
#include "nml/Group.hh"
#include "nml/Items.hh"
#include "c4/global.hh"
#include <iostream>
#include <iomanip>
#include <string>


namespace
{

using rtt_LAMGDiffusionSolver_test::multiply;

class Test_DB {
  public:

#include ".nml_test.hh"

  public:
    void setup_namelist( NML_Group& g )
    {
#include ".nml_test.cc"
    }
};

inline rtt_LAMG::Options::Outlevel outlevel(int levout)
{
    using rtt_LAMG::Options;
    
    switch (levout)
    {
    case Options::SILENT:
	return Options::SILENT;
    case Options::LEVEL1:
	return Options::LEVEL1;
    case Options::LEVEL2:
	return Options::LEVEL2;
    case Options::LEVEL3:
	return Options::LEVEL3;
    }
    Insist(false, "Illegal levout");
    return Options::SILENT;
}

}

int main( int argc, char *argv[] )
{
    std::cout << "Calling C4::Init(...)" << std::endl;
	
    C4::Init(argc, argv);
    
    using rtt_diffusion::TestSolver;
    using rtt_LAMGDiffusionSolver::SolverP1Diff;
    typedef TestSolver<Mesh_XYZ,SolverP1Diff> Tester;

    NML_Group g("testLAMGDiffusionSolver");

    rtt_diffusion::Diffusion_DB diff_db;
    diff_db.setup_namelist(g);

    Test_DB tdb;
    tdb.setup_namelist(g);
    
    Mesh_XYZ::Mesh_DB mdb;
    
    mdb.setup_namelist(g);
    
    g.readgroup("testLAMGDiffusionSolver.in");
    g.writegroup("testLAMGDiffusionSolver.out");

    rtt_dsxx::SP<Mesh_XYZ> spmesh(new Mesh_XYZ(mdb));

    const Mesh_XYZ::FieldConstructor &fCtor = spmesh;

    rtt_LAMG::Options opts;

    rtt_LAMG::Options::Outlevel levout = outlevel(tdb.levout);
    
    SolverP1Diff
	solver(opts.itsmax(tdb.itsmax).tol(tdb.tol).levout(levout));
    
    Tester tester("LAMGDiffusionSolver::SolverP1Diff", argc, argv,
		  rtt_LAMGDiffusionSolver::release(),
		  fCtor, solver, diff_db, tdb.tolerance, std::cout);
    tester.run();

    std::cout << "Calling C4::Finalize()" << std::endl;
    C4::Finalize();
    
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of testSolver.cc
//---------------------------------------------------------------------------//
