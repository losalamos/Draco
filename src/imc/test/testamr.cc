//----------------------------------*-C++-*----------------------------------//
// testamr.cc
// Thomas M. Evans
// Tue Jun  9 18:42:16 1998
//---------------------------------------------------------------------------//
// @> test parallel building and mixing for AMR IMC
//---------------------------------------------------------------------------//

#include "../AMR_Interface.hh"
#include "../AMR_Builder.hh"
#include "../OS_Mesh.hh"
#include "../Global.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <cmath> 

using IMC::AMR_Interface;
using IMC::AMR_Builder;
using IMC::OS_Mesh;
using dsxx::SP;
using namespace std;
using namespace C4;

// declare node
int mynode;
int mynodes;

int main(int argc, char *argv[])
{    
  // init C4 stuff
    C4::Init(argc, argv);
    mynode  = C4::node();
    mynodes = C4::nodes();
 
 // try block
    try
    {
      // material temp
	double tev[6]  = {10,1e-2,1e-2,10,1e-2,1e-2};
	double factor  = .20;
	double dt      = .000001;
	double max_dt  = .01;
	double elapsed = 0.0;
	for (int cycle = 0; cycle < 250; cycle++)
	{
	    
	  // let's give a mesh definition, RAGE style
	    int lay[]     = {1,1,1,2,1,4,1,2,2,3,2,5,2,3,3,3,3,6,
			     4,4,1,5,4,4,4,5,2,6,5,5,5,6,3,6,6,6};
	    double vert[] = {0,0,0,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,1, 
			     1,0,0,1,1,0,2,1,0,2,0,0,1,0,1,1,1,1,2,1,1,2,0,1,
			     2,0,0,2,1,0,3,1,0,3,0,0,2,0,1,2,1,1,3,1,1,3,0,1,
			     0,0,1,0,1,1,1,1,1,1,0,1,0,0,2,0,1,2,1,1,2,1,0,2,
			     1,0,1,1,1,1,2,1,1,2,0,1,1,0,2,1,1,2,2,1,2,2,0,2,
			     2,0,1,2,1,1,3,1,1,3,0,1,2,0,2,2,1,2,3,1,2,3,0,2};
	    for (int i = 0; i < 144; i++)
		vert[i] = vert[i] / 10.0;

	  // other data
	    int b_proc[]        = {0};
	    int b_cell[]        = {0};
	    int global_cell[]   = {1,2,3,4,5,6};
	    double dedt[6]      = {.3,.3,.3,.3,.3,.3};
	    double rho[6]       = {3,3,3,3,3,3};
	    double opacity[6]   = {100,100,100,100,100,100};
	    double rev[6]       = {10,1e-2,1e-2,10,1e-2,1e-2};
	    double t4_slope[18] = {0};
	    
	  // variables
	    int numcells    = 6;
	    int g_numcells  = 6;
	    int numbcells   = 0;
	    double implicit = 1.0;
	    double dnpdt    = 0;
	    int npnom       = 1000;
	    int npmax       = 1000;
	    int seed        = 9347593;
	    int buffer      = 1000;

	  // calculate t4-slopes
	    double f1;
	    double f2;
	    int t4i = 0;

	  // cell 1
	  // x
	    f1 = pow(tev[0], 4);
	    f2 = (pow(tev[0], 4) + pow(tev[1], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.1 - 0);
	  // y
	    t4_slope[t4i++] = 0;
	  // z
	    f2 = (pow(tev[0], 4) + pow(tev[3], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.1 - 0);

	  // cell 2
	  // x
	    f1 = (pow(tev[0], 4) + pow(tev[1], 4)) / 2;
	    f2 = (pow(tev[1], 4) + pow(tev[2], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.2 - .1);
	  // y
	    t4_slope[t4i++] = 0;
	  // z
	    f1 = pow(tev[1], 4);
	    f2 = (pow(tev[1], 4) + pow(tev[4], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.1 - 0);

	  // cell 3
	  // x
	    f1 = (pow(tev[1], 4) + pow(tev[2], 4)) / 2;
	    f2 = pow(tev[2], 4);
	    t4_slope[t4i++] = (f2 - f1) / (.3 - .2);
	  // y
	    t4_slope[t4i++] = 0;
	  // z
	    f1 = pow(tev[2], 4);
	    f2 = (pow(tev[2], 4) + pow(tev[5], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.1 - 0);

	  // cell 4
	  // x
	    f1 = pow(tev[3], 4);
	    f2 = (pow(tev[3], 4) + pow(tev[4], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.1 - 0);
	  // y
	    t4_slope[t4i++] = 0;
	  // z
	    f1 = (pow(tev[0], 4) + pow(tev[3], 4)) / 2;
	    f2 = pow(tev[3], 4);
	    t4_slope[t4i++] = (f2 - f1) / (.2 - .1);

	  // cell 5
	  // x
	    f1 = (pow(tev[3], 4) + pow(tev[4], 4)) / 2;
	    f2 = (pow(tev[4], 4) + pow(tev[5], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.2 - .1);
	  // y
	    t4_slope[t4i++] = 0;
	  // z
	    f1 = (pow(tev[1], 4) + pow(tev[4], 4)) / 2;
	    f2 = pow(tev[4], 4);
	    t4_slope[t4i++] = (f2 - f1) / (.2 - .1);

	  // cell 6
	  // x
	    f2 = pow(tev[5], 4);
	    f1 = (pow(tev[4], 4) + pow(tev[5], 4)) / 2;
	    t4_slope[t4i++] = (f2 - f1) / (.3 - .2);
	  // y
	    t4_slope[t4i++] = 0;
	  // z
	    f1 = (pow(tev[2], 4) + pow(tev[5], 4)) / 2;
	    f2 = pow(tev[5], 4);
	    t4_slope[t4i++] = (f2 - f1) / (.2 - .1);
	    Check (t4i == 18);
	    
	  // return arrays
	    double e_dep[6]   = {0};
	    double rad_den[6] = {0};

	  // initialize
	    for (int i = 0; i < numcells; i++)
		opacity[i] = opacity[i] / (tev[i] * tev[i] * tev[i]);

	  // call rage_imc
	    rage_imc_(&cycle, &numcells, &g_numcells, vert, lay, &numbcells,
		      b_proc, b_cell, dedt, tev, rev, rho, opacity, opacity,
		      &dt, &elapsed, &implicit, &npnom, &npmax, &dnpdt,
		      &seed, &buffer, global_cell, t4_slope, e_dep, rad_den);
	    
	  // update temp
	    cout << " *** RESULTS FOR CYCLE " << cycle << " ***" << endl;
	    cout << "  ** Elapsed Time : " << elapsed << endl;
	    cout << "  ** Time Step    : " << dt << endl;
	    for (int i = 0; i < numcells; i++)
	    {
		tev[i] += e_dep[i] / dedt[i];
		cout << i+1 << "\t" << tev[i] << "\t" << rad_den[i] << endl;
	    }
	    cout << endl;
	    
	  // increment timestep
	    if (cycle > 0 && dt < max_dt)
		dt += dt * factor;
	    elapsed += dt;
	}
    }
    catch (const dsxx::assertion &ass)
    {
	cout << "Dumbass, you screwed up: " << ass.what() << endl;
	return 1;
    }
    catch(...)
    {
	cout << "HELP ME" << endl;
	return 1;
    }

  // c4 end
    C4::Finalize();

  // we ran ok
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of testamr.cc
//---------------------------------------------------------------------------//
