//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Fri Mar 20 12:04:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/testFullP13T.hh"
#include "mesh/Mesh_XYZ.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"

#include <stdexcept>
#include <iostream>
#include <string>

using std::cerr;
using std::cout;
using std::endl;

typedef Mesh_XYZ MT;

namespace {

 const std::string infile = "testFullP13T.in";
 const std::string outfile = "testFullP13T.out";

 template<class UMCMP>
 void runProblem(const XTM::testFullP13T_DB &tdb, const Diffusion_DB &diffdb,
		 const Mesh_DB &mdb, const pcg_DB &pcg_db)
 {
     XTM::testFullP13T<MT,UMCMP> test(tdb, diffdb, mdb, pcg_db);

     cerr << C4::node() << " before test.run()" << endl;

     test.run();

     cerr << C4::node() << " after test.run()" << endl;
 }

 int run(int argc, char *argv[])
 {
     try
     {
	 NML_Group g("testFullP13T");

	 XTM::testFullP13T_DB tdb;
	 tdb.setup_namelist(g);

	 Diffusion_DB diffdb;
	 diffdb.setup_namelist(g);

	 Mesh_DB mdb;
	 mdb.setup_namelist(g);

	 pcg_DB pcg_db("pcg");
	 pcg_db.setup_namelist(g);

	 if (C4::node() == 0)
	     cout << "Reading input from " << infile << ".\n";

	 g.readgroup(infile.c_str());
	 g.writegroup(outfile.c_str());

	 if (C4::node() == 0)
	     cerr << "In main()" << endl;
	 
	 switch (tdb.matprops)
	 {
	 case XTM::Marshak:
	     runProblem<rtt_matprops::MarshakMaterialProps>(tdb, diffdb, mdb,
							    pcg_db);
	     break;
	 case XTM::Interped:
	     runProblem<rtt_matprops::InterpedMaterialProps>(tdb, diffdb, mdb,
							     pcg_db);
	     break;
	 default:
	     runProblem<rtt_matprops::MarshakMaterialProps>(tdb, diffdb, mdb,
							    pcg_db);
	     break;
	 }

     }
     catch (const char *str)
     {
	 cerr << "caught: " << str << endl;
	 return 1;
     }
     catch (const dsxx::assertion &ass)
     {
	 cerr << "caught assertion exception: " << ass.what() << endl;
	 return 1;
     }
     catch (const std::runtime_error &rtx)
     {
	 cerr << "caught assertion exception: " << rtx.what() << endl;
	 return 1;
     }
     catch (...)
     {
	 cerr << "caught unknown exception" << endl;
	 return 1;
     }
     return 0;
 }

} // end unnamed namespace

int main(int argc, char *argv[])
{
    C4::Init( argc, argv );

    {
	// Introduce a new scope in case the test problem object may need
	// to do MPI work in its dtor.

	run(argc, argv);
    }

    C4::Finalize();    
}

#include "3T/testP13T/testFullP13T.cc"

template
class XTM::testFullP13T<MT, rtt_matprops::InterpedMaterialProps>;


template
class XTM::testFullP13T<MT, rtt_matprops::MarshakMaterialProps>;

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
