//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Thu Sep 24 16:45:25 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "testP1Diffusion.hh"

#include "../DiffusionSelector.hh"
#include "nml/Group.hh"
#include "nml/Item.hh"
#include "nml/Items.hh"
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

	typedef rtt_P1Diffusion::DiffusionSelector<MT>::Options Options;
	
#ifdef SELECTOR_PCG
	Options options("pcg");
#endif
	mdb.setup_namelist(g);
    
#ifdef SELECTOR_PCG
	options.setup_namelist(g);
#endif
	g.readgroup("testP1Diffusion.in");
	g.writegroup("testP1Diffusion.out");

	mdb.resize();
	
#ifdef SELECTOR_CONJGRAD
	Options options(tdb.maxIterations, tdb.epsilon);
#endif	

	SP<MT> spMesh(new MT(mdb));
	rtt_P1Diffusion_test::testP1Diffusion<MT> testP1(spMesh, spMesh,
							 tdb.D, tdb.sigma,
							 tdb.q, tdb.fTop,
							 tdb.fBot,
							 diffdb, options);
	double error = testP1.run();

	if (error <= tdb.tolerance)
	    std::cout << "**** " << "testP1Diffusion"
		      << " Test: PASSED ****" << std::endl;
	else
	    std::cout << "**** The solver error, " << error
		      << ", did not meet tolerance, " << tdb.tolerance << "."
		      << " ****" << std::endl
		      << "**** " << "testP1Diffusion"
		      << " Test: FAILED ****" << std::endl;
 
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
