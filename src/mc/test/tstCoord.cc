//----------------------------------*-C++-*----------------------------------//
// tstCoord.cc
// Thomas M. Evans
// Fri Apr 16 14:33:48 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Test of MC Coord_sys classes 
//---------------------------------------------------------------------------//

#include "../Coord_sys.hh"
#include "../XYCoord_sys.hh"
#include "../XYZCoord_sys.hh"
#include "../Release.hh"
#include "rng/Rnd_Control.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_mc::Coord_sys;
using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys;

using dsxx::SP;

using namespace std;

// some necessary global stuff
vector<double> random_nums(100, 0.0);
int seed = 9375632;
Rnd_Control rndc(seed);
bool passed = true;

// make some random numbers for testing
void make_ran()
{
    // make a Sprng
    Sprng x = rndc.get_rn(10);

    // fill up the array
    for (int i = 0; i < 100; i++)
	random_nums[i] = x.ran();
}

// test XYCoord
void testXY(Coord_sys &xy)
{
    // test dimension functions
    if (xy.get_dim() != 2)      passed = false;
    if (xy.get_sdim() != 3)     passed = false;
    if (xy.get_Coord() != "xy") passed = false;

    // test virtual functions
    
    // test a sampling of position in a cell
    Sprng ts(rndc.get_rn(10));
    vector<double> min(2, -1);
    vector<double> max(2, 1);
    vector<double> answer = xy.sample_pos(min, max, ts);
    vector<double> control(2, 0.0);
    control[0] = 2 * random_nums[0] + -1;
    control[1] = 2 * random_nums[1] + -1;
    if (answer != control) passed = false;
}

// test XYZCoord
void testXYZ(Coord_sys &xyz)
{
    // test dimension functions
    if (xyz.get_dim() != 3)       passed = false;
    if (xyz.get_sdim() != 3)      passed = false;
    if (xyz.get_Coord() != "xyz") passed = false;
}

int main(int argc, char *argv[])
{  
    C4::Init(argc, argv);

    // this is a scalar test
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

    // fill up the random number array
    make_ran();

    // make some Coord_sys stuff
    XYCoord_sys xy;
    XYZCoord_sys xyz;
    SP<XYCoord_sys> sxy(new XYCoord_sys());
    SP<XYZCoord_sys> sxyz(new XYZCoord_sys());

    // test XY
    testXY(xy);
    testXY(*sxy);

    // test XYZ
    testXYZ(xyz);
    testXYZ(*sxyz);

    // status of test
    cout << endl;
    cout <<     "*************************************" << endl;
    if (passed) 
    {
        cout << "**** Coord_sys Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** Coord_sys Self Test: FAILED ****" << endl;
    }
    cout <<     "*************************************" << endl;
    cout << endl;

    cout << "Done testing Coord_sys." << endl;

    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                              end of tstCoord.cc
//---------------------------------------------------------------------------//
