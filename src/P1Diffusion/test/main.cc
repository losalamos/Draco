//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Thu Sep 24 16:45:25 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "testP1Diffusion.hh"

#include "PCGDiffusionSolver/pcg_DB.hh"
#include "nml/Group.hh"
#include "ds++/SP.hh"
#include "c4/global.hh"
#include "mesh/Mesh_XYZ.hh"
#include <iostream>
#include <iomanip>
#include <string>

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

using rtt_dsxx::SP;

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    std::cout << progname << ": version " << version << std::endl;
}

int main(int argc, char *argv[])
{    
    C4::Init( argc, argv );

    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    C4::Finalize();
	    return 0;
	}
    }

    try
    {

	std::cout << "Initiating " << "testP1Diffsuion" << " test"
		  << std::endl;

	NML_Group g("testP1Diffusion");

	rtt_diffusion::Diffusion_DB diffdb;
	diffdb.setup_namelist(g);

	Test_DB tdb;
	tdb.setup_namelist(g);
    
	rtt_mesh::Mesh_DB mdb;

	using rtt_PCGDiffusionSolver::pcg_DB;
	pcg_DB pcg_db("pcg");
    
	mdb.setup_namelist(g);
    
	pcg_db.setup_namelist(g);

	g.readgroup("testP1Diffusion.in");
	g.writegroup("testP1Diffusion.out");

	SP<MT> spMesh(new MT(mdb));
	rtt_P1Diffusion_test::testP1Diffusion<MT> testP1(spMesh, spMesh,
							 tdb.D, tdb.sigma,
							 tdb.q, tdb.fTop,
							 tdb.fBot,
							 diffdb, pcg_db);
	testP1.run();

	std::cout << "**** " << "testP1Diffusion"
		  << " Test: PASSED ****" << std::endl;
 
    }
    catch( rtt_dsxx::assertion& a )
    {
	std::cout << "Failed assertion: " << a.what() << std::endl;
	std::cout << "**** " << "testP1Diffusion"
		  << " Test: FAILED ****" << std::endl;
    }
    catch( std::exception& a )
    {
	std::cout << "Failed exception: " << a.what() << std::endl;
	std::cout << "**** " << "testP1Diffusion"
		  << " Test: FAILED ****" << std::endl;
    }
    
    C4::Finalize();

    return 0;
}

#include "testP1Diffusion.t.hh"
template
class rtt_P1Diffusion_test::testP1Diffusion<MT>;

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
