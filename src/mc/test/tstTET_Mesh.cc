//----------------------------------*-C++-*----------------------------------//
// tstTET_Mesh.cc
// H. Grady Hughes
// Thu Jan 27 16:34:07 MST 2000
// $Id$
//---------------------------------------------------------------------------//
// @> Test of TET_Mesh and TET_Builder
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../TET_Mesh.hh"
#include "../Layout.hh"
#include "../XYCoord_sys.hh"
#include "../XYZCoord_sys.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_mc::Layout;
using rtt_mc::TET_Mesh;
using dsxx::SP;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

// mesh proxy class
class Mesh_Proxy
{
  private:
    SP<TET_Mesh> mesh;
  public:
    Mesh_Proxy(SP<TET_Mesh> m) : mesh(m) {}
    const TET_Mesh& get_Mesh() const { return *mesh; }
};

// 2D Mesh Tests
void Test_2D()
{
}

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_mc::release() << endl;
	    C4::Finalize();
	    return 0;
	}

    // 2D Mesh tests
    Test_2D();

    // status of test
    cout << endl;
    cout <<     "************************************" << endl;
    if (passed) 
    {
        cout << "**** TET_Mesh Self Test: PASSED ****" << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing TET_Mesh." << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                              end of tstTET_Mesh.cc
//---------------------------------------------------------------------------//
