//----------------------------------*-C++-*----------------------------------//
// testamr.cc
// Thomas M. Evans
// Tue Jun  9 18:42:16 1998
//---------------------------------------------------------------------------//
// @> test parallel building and mixing for AMR IMC
//---------------------------------------------------------------------------//

#include "imc/AMR_Interface.hh"
#include "imc/AMR_Builder.hh"
#include "imc/OS_Mesh.hh"
#include "imc/Global.hh"
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
	double tev[6] = {1,1,1,1,1,1};
	for (int cycle = 0; cycle < 10; cycle++)
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

	  // other data
	    int b_proc[]       = {0};
	    int b_cell[]       = {0};
	    int global_cell[]  = {1,2,3,4,5,6};
	    double dedt[6]     = {.3,.3,.3,.3,.3,.3};
	    double rho[6]      = {3,3,3,3,3,3};
	    double opacity[6]  = {100,100,100,100,100,100};
	    double rev[6]      = {1,1,1,1,1,1};
	    double t4_slope[6] = {0};
	    
	  // variables
	    int numcells    = 6;
	    int g_numcells  = 6;
	    int numbcells   = 0;
	    double implicit = 1.0;
	    double dt       = .001;
	    double dnpdt    = 0;
	    int npnom       = 1000;
	    int npmax       = 1000;
	    int seed        = 9347593;
	    int buffer      = 1000;
	    
	  // return arrays
	    double e_dep[6]   = {0};
	    double rad_den[6] = {0};

	  // initialize
	    for (int i = 0; i < numcells; i++)
		opacity[i] = opacity[i] / (tev[i] * tev[i] * tev[i]);
	    double elapsed = IMC::Global::max(0.0, (cycle - 1) * dt);

	  // call rage_imc
	    rage_imc_(&cycle, &numcells, &g_numcells, vert, lay, &numbcells,
		      b_proc, b_cell, dedt, tev, rev, rho, opacity, opacity,
		      &dt, &elapsed, &implicit, &npnom, &npmax, &dnpdt,
		      &seed, &buffer, global_cell, t4_slope, e_dep, rad_den);
	    
	  // update temp
	    cout << " *** RESULTS FOR CYCLE " << cycle << " ***" << endl;
	    for (int i = 0; i < numcells; i++)
	    {
		tev[i] += e_dep[i] / dedt[i];
		cout << i+1 << "\t" << tev[i] << "\t" << rad_den[i] << endl;
	    }
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
