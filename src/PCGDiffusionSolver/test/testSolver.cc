//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCGDiffusionSolver/test/testSolver.cc
 * \author Randy M. Roberts
 * \date   Mon Feb  7 10:49:55 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../SolverP1Diff.t.hh"
#include "../MatrixP1Diff.t.hh"
#include "../MatVecP1Diff.t.hh"
#include "../PreCondP1Diff.t.hh"
#include "../MatrixP1DiffTraits.hh"
#include "../Release.hh"
#include "../pcg_DB.hh"

#include "diffusion/TestSolver.t.hh"
#include "mesh/Mesh_XYZ.hh"
#include "nml/Group.hh"
#include <iostream>
#include <iomanip>
#include <string>

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

}

int main( int argc, char *argv[] )
{
    using rtt_diffusion::TestSolver;
    using rtt_PCGDiffusionSolver::SolverP1Diff;
    using rtt_PCGDiffusionSolver::pcg_DB;

    typedef TestSolver<Mesh_XYZ,SolverP1Diff<Mesh_XYZ> > Tester;

    NML_Group g("testPCGDiffusionSolver");

    rtt_diffusion::Diffusion_DB diff_db;
    diff_db.setup_namelist(g);

    Test_DB tdb;
    tdb.setup_namelist(g);
    
    Mesh_DB mdb;
    pcg_DB pcg_db("pcg");
    
    mdb.setup_namelist(g);
    
    pcg_db.setup_namelist(g);

    g.readgroup("testPCGDiffusionSolver.in");
    g.writegroup("testPCGDiffusionSolver.out");

    rtt_dsxx::SP<Mesh_XYZ> spmesh(new Mesh_XYZ(mdb));

    const Mesh_XYZ::FieldConstructor &fCtor = spmesh;

    SolverP1Diff<Mesh_XYZ> solver(spmesh, pcg_db);
    
    Tester tester("PCGDiffusionSolver::SolverP1Diff", argc, argv,
		  rtt_PCGDiffusionSolver::release(),
		  fCtor, solver, diff_db, tdb.tolerance, std::cout);
    tester.run();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of testSolver.cc
//---------------------------------------------------------------------------//
