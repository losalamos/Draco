//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstCoord.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 20 16:28:17 2001
 * \brief  Coord_sys test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../Coord_sys.hh"
#include "../XYCoord_sys.hh"
#include "../XYZCoord_sys.hh"
#include "../Math.hh"
#include "../Release.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_mc::Coord_sys;
using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys; 
using rtt_mc::global::dot;

using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

// some necessary global stuff
vector<double> random_nums(100, 0.0);
int seed = 9375632;
Rnd_Control rndc(seed);


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// make some random numbers for testing
void make_ran()
{
    // make a Sprng
    Sprng x = rndc.get_rn(10);

    // fill up the array
    for (int i = 0; i < 100; i++)
	random_nums[i] = x.ran();
}

//---------------------------------------------------------------------------//

// test XYCoord
void testXY(Coord_sys &xy)
{
    // test dimension functions
    if (xy.get_dim() != 2)      ITFAILS;
    if (xy.get_sdim() != 3)     ITFAILS;
    if (xy.get_Coord() != "xy") ITFAILS;

    // test virtual functions
    
    // test a sampling of position in a cell
    Sprng ts(rndc.get_rn(10));
    vector<double> min(2, -1.0);
    vector<double> max(2, 1.0);
    vector<double> answer = xy.sample_pos(min, max, ts);
    vector<double> control(2, 0.0);
    control[0] = 2 * random_nums[0] + -1;
    control[1] = 2 * random_nums[1] + -1;
    
    for (int i = 0; i < answer.size(); i++)
	if (!soft_equiv(answer[i], control[i])) ITFAILS;
}

//---------------------------------------------------------------------------//

// test XYZCoord
void testXYZ(Coord_sys &xyz)
{
    // test dimension functions
    if (xyz.get_dim() != 3)       ITFAILS;
    if (xyz.get_sdim() != 3)      ITFAILS;
    if (xyz.get_Coord() != "xyz") ITFAILS;
}

//---------------------------------------------------------------------------//

// test the direction cosines from initial sampling through several
// scatters. 

void test_scatters(Coord_sys &xyz)
{
    // direction vector
    vector<double> direction(3, 0.0);

    // 
    // test scatter back and forth along the x-axis
    //

    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;

    for (int i = 0; i < 1000; i++)
    {
	xyz.calc_omega(-1.0, 0.0, direction);
	
	if (!soft_equiv(direction[0], pow(-1.0,i+1), 1.e-10)) ITFAILS;
	if (!soft_equiv(direction[1],           0.0, 1.e-10)) ITFAILS;
	if (!soft_equiv(direction[2],           0.0, 1.e-10)) ITFAILS;

	if (!soft_equiv(sqrt(dot(direction,direction)),1.0,1.e-10)) ITFAILS; 
    }

    // 
    // test scatter back and forth along the y-axis
    //

    direction[0] = 0.0;
    direction[1] = 1.0;
    direction[2] = 0.0;

    for (int i = 0; i < 1000; i++)
    {
	xyz.calc_omega(-1.0, 0.0, direction);
	
	if (!soft_equiv(direction[0],           0.0, 1.e-10)) ITFAILS;
	if (!soft_equiv(direction[1], pow(-1.0,i+1), 1.e-10)) ITFAILS;
	if (!soft_equiv(direction[2],           0.0, 1.e-10)) ITFAILS;

	if (!soft_equiv(sqrt(dot(direction,direction)),1.0,1.e-10)) ITFAILS; 
    }

    // 
    // test scatter back and forth along the z-axis
    //

    direction[0] = 0.0;
    direction[1] = 0.0;
    direction[2] = 1.0;

    for (int i = 0; i < 1000; i++)
    {
	xyz.calc_omega(-1.0, 0.0, direction);
	
	if (!soft_equiv(direction[0],           0.0, 1.e-10)) ITFAILS;
	if (!soft_equiv(direction[1],           0.0, 1.e-10)) ITFAILS;
	if (!soft_equiv(direction[2], pow(-1.0,i+1), 1.e-10)) ITFAILS;

	if (!soft_equiv(sqrt(dot(direction,direction)),1.0,1.e-10)) ITFAILS; 
    }

    //
    // begin with an isotropic direction and scatter many times (robustness) 
    //

    // make a Sprng random number object
    Sprng x = rndc.get_rn(10);

    // get an initial direction
    direction = xyz.sample_dir("isotropic", x);

    // make several scatters through a constant angle
    for (int i = 0; i < 1000; i++)
	xyz.calc_omega(0.5,0.75,direction);

    //
    // begin anew along the z-axis, and scatter many times (robustness)
    //

    direction[0] = 0.0;
    direction[1] = 0.0;
    direction[2] = 1.0;

    // make several scatters through a constant angle
    for (int i = 0; i < 1000; i++)
	xyz.calc_omega(0.5,0.75,direction);
}

//---------------------------------------------------------------------------//

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
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
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
	
	// test direction cosines with scattering
	test_scatters(xyz);
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstCoord, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstCoord Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstCoord on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstCoord.cc
//---------------------------------------------------------------------------//
