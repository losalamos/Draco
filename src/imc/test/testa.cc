//----------------------------------*-C++-*----------------------------------//
// testa.cc
// Thomas M. Evans
// Tue Jun  9 18:42:16 1998
//---------------------------------------------------------------------------//
// @> test parallel building and mixing for AMR IMC
//---------------------------------------------------------------------------//

#include "imc/AMR_Interface.hh"
#include "imc/AMR_Builder.hh"
#include "imc/Opacity_Builder.hh"
#include "imc/Mat_State.hh"
#include "imc/Opacity.hh"
#include "imc/Parallel_Source_Init.hh"
#include "imc/OS_Mesh.hh"
#include "imc/Global.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
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
using IMC::Opacity_Builder;
using IMC::Mat_State;
using IMC::Opacity;
using IMC::Parallel_Source_Init;
using RNG::Rnd_Control;
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
      // declare geometry and material stuff
	SP<OS_Mesh> mesh;
	SP<Mat_State<OS_Mesh> > mat_state;
	SP<Opacity<OS_Mesh> > opacity;
	SP<Parallel_Source_Init<OS_Mesh> > source_init;
	SP<Rnd_Control> rnd_con = new Rnd_Control(3423432);

	for (int cycle = 0; cycle < 2; cycle++)
	{
	    
	  // let's give a mesh definition, RAGE style
	    int lay[]     = {0,0,0,2,0,0,1,0,0,0,0,0};
	    double vert[] = {0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1,
			     1,0,0,2,0,0,2,1,0,1,1,0,1,0,1,2,0,1,2,1,1,1,1,1};
	    int b1[]      = {0};
	    int b2[]      = {0};
	    int num_cells = 2;
	    int num_b     = 0;
	    int gc[]      = {1,2};
	    double cv[]   = {.1,.1};
	    double real   = 1;
	    int integer   = 10;
	    
	  // let's get an interface to this problem
	    cout << ">> Running through problem builders" << endl;
	    
	    AMR_Interface::Arguments arg(vert, lay, b1, b2, gc, cv, cv, cv, cv,
					 cv, num_cells, num_cells, num_b, real,
					 real, real, real, integer, integer,
					 real, integer, integer, integer, 
					 cycle);
	    SP<AMR_Interface> interface = new AMR_Interface(arg);
	    cout << "** Merged with host data on node " << mynode << endl;
	    
	  // build the mesh
	    AMR_Builder amr_build(interface);
	    mesh = amr_build.build_Mesh();
	    cout << "** Built mesh on node " << mynode << endl;
	    
	  // build the Opacities
	    Opacity_Builder<OS_Mesh> opacity_builder(interface);
	    mat_state = opacity_builder.build_Mat(mesh);
	    opacity   = opacity_builder.build_Opacity(mesh, mat_state);
	    cout << "** Built opacities on node " << mynode << endl;

	  // do the source initialization
	    source_init = new Parallel_Source_Init<OS_Mesh>(interface, mesh);
	    if (cycle == 0)
		AMR_Interface::set_census(source_init->calc_initial_census
					  (mesh, opacity, mat_state,
					   rnd_con));
	    cout << "** We did source initialization on node " << mynode
		 << endl;
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
//                              end of testa.cc
//---------------------------------------------------------------------------//
