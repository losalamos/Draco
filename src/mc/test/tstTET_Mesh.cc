//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstTET_Mesh.cc
 * \author H. Grady Hughes
 * \date   Thu Jan 27 16:34:07 MST 2000
 * \brief  Tests for tetrahedral meshes and ThreeVectors.
 */
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../TET_Mesh.hh"
#include "../Layout.hh"
#include "../XYCoord_sys.hh"
#include "../XYZCoord_sys.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"

#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_mc::Layout;
using rtt_mc::TET_Mesh;
using rtt_mc::ThreeVector;
using rtt_mc::sample_in_triangle;
using rtt_rng::Sprng;
using dsxx::SP;

int seed = 493875348;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

//! Mesh proxy class
class Mesh_Proxy
{
  private:
    SP<TET_Mesh> mesh;
  public:
    Mesh_Proxy(SP<TET_Mesh> m) : mesh(m) {}
    const TET_Mesh& get_Mesh() const { return *mesh; }
};

// ThreeVector Tests
void Test_ThreeVector()
{
    int num = 5;
    int *idr = init_sprng(0, num, seed, 1);
    Sprng random(idr, 0);

    // Test sampling in a triangle.

    int numm = 10000;
    double factor = 36.0/static_cast<double>(numm);
    double eps_triang = 0.0001;

    ThreeVector A(100.0,10.0,0.0);
    ThreeVector B(112.0,10.0,0.0);
    ThreeVector C(106.0,16.0,0.0);

    double count[6][12];
    for (int y = 0; y < 6 ; y++)
        for (int x = 0; x < 12 ; x++)
            count[y][x] = 0.;

    for (int i = 0 ; i < numm ; i++)
    {
        ThreeVector R = sample_in_triangle(A,B,C,random);

        if (R.get_x() <= 100.0 || R.get_x() >= 112.0)       ITFAILS;
        if (R.get_y() <= 10.0  || R.get_y() >= 16.0)        ITFAILS;
        if (R.get_z() != 0.0)                               ITFAILS;

        int xx = static_cast<int>(R.get_x() - 100.);
        int yy = static_cast<int>(R.get_y() - 10.);
        count[yy][xx] += 1.;
    }

    for (int yy = 5 ; yy >= 0 ; yy--)
        for (int xx = 0; xx < 12; xx++)
            count[yy][xx] *= factor;

    if (count[5][0] != 0.0)                        ITFAILS;
    if (count[5][1] != 0.0)                        ITFAILS;
    if (count[5][2] != 0.0)                        ITFAILS;
    if (count[5][3] != 0.0)                        ITFAILS;
    if (count[5][4] != 0.0)                        ITFAILS;
    if (fabs(count[5][5] - 0.5544) > eps_triang)   ITFAILS;
    if (fabs(count[5][6] - 0.5040) > eps_triang)   ITFAILS;
    if (count[5][7] != 0.0)                        ITFAILS;
    if (count[5][8] != 0.0)                        ITFAILS;
    if (count[5][9] != 0.0)                        ITFAILS;
    if (count[5][10] != 0.0)                       ITFAILS;
    if (count[5][11] != 0.0)                       ITFAILS;
    if (count[4][0] != 0.0)                        ITFAILS;
    if (count[4][1] != 0.0)                        ITFAILS;
    if (count[4][2] != 0.0)                        ITFAILS;
    if (count[4][3] != 0.0)                        ITFAILS;
    if (fabs(count[4][4] - 0.4968) > eps_triang)   ITFAILS;
    if (fabs(count[4][5] - 1.0476) > eps_triang)   ITFAILS;
    if (fabs(count[4][6] - 1.0368) > eps_triang)   ITFAILS;
    if (fabs(count[4][7] - 0.5040) > eps_triang)   ITFAILS;
    if (count[4][8] != 0.0)                        ITFAILS;
    if (count[4][9] != 0.0)                        ITFAILS;
    if (count[4][10] != 0.0)                       ITFAILS;
    if (count[4][11] != 0.0)                       ITFAILS;
    if (count[3][0] != 0.0)                        ITFAILS;
    if (count[3][1] != 0.0)                        ITFAILS;
    if (count[3][2] != 0.0)                        ITFAILS;
    if (fabs(count[3][3] - 0.5220) > eps_triang)   ITFAILS;
    if (fabs(count[3][4] - 0.9792) > eps_triang)   ITFAILS;
    if (fabs(count[3][5] - 0.9396) > eps_triang)   ITFAILS;
    if (fabs(count[3][6] - 0.9504) > eps_triang)   ITFAILS;
    if (fabs(count[3][7] - 1.0260) > eps_triang)   ITFAILS;
    if (fabs(count[3][8] - 0.5184) > eps_triang)   ITFAILS;
    if (count[3][9] != 0.0)                        ITFAILS;
    if (count[3][10] != 0.0)                       ITFAILS;
    if (count[3][11] != 0.0)                       ITFAILS;
    if (count[2][0] != 0.0)                        ITFAILS;
    if (count[2][1] != 0.0)                        ITFAILS;
    if (fabs(count[2][2] - 0.4932) > eps_triang)   ITFAILS;
    if (fabs(count[2][3] - 0.9972) > eps_triang)   ITFAILS;
    if (fabs(count[2][4] - 1.0800) > eps_triang)   ITFAILS;
    if (fabs(count[2][5] - 1.0872) > eps_triang)   ITFAILS;
    if (fabs(count[2][6] - 1.0404) > eps_triang)   ITFAILS;
    if (fabs(count[2][7] - 0.9324) > eps_triang)   ITFAILS;
    if (fabs(count[2][8] - 1.0368) > eps_triang)   ITFAILS;
    if (fabs(count[2][9] - 0.5508) > eps_triang)   ITFAILS;
    if (count[2][10] != 0.0)                       ITFAILS;
    if (count[2][11] != 0.0)                       ITFAILS;
    if (count[1][0] != 0.0)                        ITFAILS;
    if (fabs(count[1][1] - 0.4752) > eps_triang)   ITFAILS;
    if (fabs(count[1][2] - 1.0116) > eps_triang)   ITFAILS;
    if (fabs(count[1][3] - 1.0044) > eps_triang)   ITFAILS;
    if (fabs(count[1][4] - 0.9684) > eps_triang)   ITFAILS;
    if (fabs(count[1][5] - 1.0188) > eps_triang)   ITFAILS;
    if (fabs(count[1][6] - 0.9324) > eps_triang)   ITFAILS;
    if (fabs(count[1][7] - 1.0080) > eps_triang)   ITFAILS;
    if (fabs(count[1][8] - 1.0044) > eps_triang)   ITFAILS;
    if (fabs(count[1][9] - 0.8424) > eps_triang)   ITFAILS;
    if (fabs(count[1][10] - 0.5076) > eps_triang)  ITFAILS;
    if (count[1][11] != 0.0)                       ITFAILS;
    if (fabs(count[0][0] - 0.5184) > eps_triang)   ITFAILS;
    if (fabs(count[0][1] - 0.9216) > eps_triang)   ITFAILS;
    if (fabs(count[0][2] - 0.9180) > eps_triang)   ITFAILS;
    if (fabs(count[0][3] - 1.0872) > eps_triang)   ITFAILS;
    if (fabs(count[0][4] - 1.0080) > eps_triang)   ITFAILS;
    if (fabs(count[0][5] - 0.9900) > eps_triang)   ITFAILS;
    if (fabs(count[0][6] - 0.9540) > eps_triang)   ITFAILS;
    if (fabs(count[0][7] - 1.0188) > eps_triang)   ITFAILS;
    if (fabs(count[0][8] - 1.0044) > eps_triang)   ITFAILS;
    if (fabs(count[0][9] - 1.0548) > eps_triang)   ITFAILS;
    if (fabs(count[0][10] - 0.9936) > eps_triang)  ITFAILS;
    if (fabs(count[0][11] - 0.4608) > eps_triang)  ITFAILS;

//  Beginning:  a way to see these results ordered on the page.
//  cout.setf(ios::fixed);
//  cout.precision(5);
//
//  for (int yy = 5 ; yy >= 0 ; yy--) {
//      for (int xx = 0; xx < 12; xx++)
//          cout << setw(10) << factor*count[yy][xx];
//      cout << endl;
//  }
//  End:  a way to see these results ordered on the page.

    // End of test sampling in a triangle.

}

// TET_Mesh Tests
void Test_TET()
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

    // ThreeVector tests
    Test_ThreeVector();

    // TET_Mesh tests
    Test_TET();

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
