//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Thu Sep 24 16:45:25 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/test/testSolverP1Diff.hh"
#include "P1Diffusion/test/testP1Diffusion.hh"

#include "nml/Group.hh"
#include "ds++/SP.hh"
#include "mesh/Mesh_XYZ.hh"
#include <iostream>
#include <iomanip>

class Test_DB {
  public:

#include ".nml_test.hh"

  public:
    void setup_namelist( NML_Group& g )
    {
#include ".nml_test.cc"
    }
};

typedef Mesh_XYZ MT;
typedef MT::ccsf ccsf;

using dsxx::SP;

int main()
{    
    NML_Group g("testP1Diffusion");

    Diffusion_DB diffdb;
    diffdb.setup_namelist(g);

    Test_DB tdb;
    tdb.setup_namelist(g);
    
    Mesh_DB mdb;
    pcg_DB pcg_db("pcg");
    
    mdb.setup_namelist(g);
    
    pcg_db.setup_namelist(g);

    g.readgroup("testP1Diffusion.in");
    g.writegroup("testP1Diffusion.out");

    rtt_P1Diffusion_test::testSolverP1Diff<MT> testSolve(mdb, pcg_db);
    testSolve.run();

    SP<MT> spMesh = new MT(mdb);
    rtt_P1Diffusion_test::testP1Diffusion<MT> testP1(spMesh, tdb.D, tdb.sigma,
						     tdb.q, tdb.fTop, tdb.fBot,
						     diffdb, pcg_db);
    testP1.run();
    
    return 0;
}

#include "testSolverP1Diff.t.cc"

template
class rtt_P1Diffusion_test::testSolverP1Diff<MT>;

#include "testP1Diffusion.t.cc"
template
class rtt_P1Diffusion_test::testP1Diffusion<MT>;

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
