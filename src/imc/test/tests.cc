//----------------------------------*-C++-*----------------------------------//
// tests.cc
// Thomas M. Evans
// Mon Jun 22 20:58:26 1998
//---------------------------------------------------------------------------//
// @> test of source initialization
//---------------------------------------------------------------------------//

#include "imc/OS_Interface.hh"
#include "imc/OS_Builder.hh"
#include "imc/OS_Mesh.hh"
#include "imc/Mat_State.hh"
#include "imc/Opacity_Builder.hh"
#include "imc/Opacity.hh"
#include "imc/Source_Init.hh"
#include "imc/Global.hh"
#include "rng/Rnd_Control.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <iomanip>
#include <string>

using namespace IMC;
using namespace std;
using namespace dsxx;

template<class MT>
void print_sinit(const MT &mesh, const Source_Init<MT> &sinit)
{
    cout.precision(3);
    cout.setf(ios::scientific);

  // print out volume emission
    cout.setf(ios::right);
    cout << setw(10) << "Cell" << setw(15) << "nvol" 
	 << setw(15) << "evol" << setw(15) << "evol_net"
	 << setw(15) << "ew_vol" << endl;
    for (int i = 1; i <= mesh.num_cells(); i++)
	cout << setw(10) << i << setw(15) << sinit.nvol(i) 
	     << setw(15) << sinit.evol(i) << setw(15) << sinit.evol_net(i) 
	     << setw(15) << sinit.ew_vol(i) << endl;
    cout << "Eloss: " << sinit.eloss_vol << endl;

    double esource = 0;
    for (int i = 1; i <= mesh.num_cells(); i++)
	esource += sinit.evol(i) + sinit.ess(i);
    cout << "Energy balance: " << setw(10) << sinit.evoltot + sinit.esstot
	 << setw(10) << esource << endl;
}

int main(int argc, char *argv[])
{
  // try block
    try
    {
      // declare geometry and material stuff
	SP<OS_Mesh> mesh;
	SP< Mat_State<OS_Mesh> > mat_state;
	SP< Opacity<OS_Mesh> > opacity;
	SP< Source_Init<OS_Mesh> > sinit;
	SP<Rnd_Control> rcon = new Rnd_Control(9836592);

      // scoping blocks for build-stuff
	{
	    string infile = argv[1];

	  // run the interface parser
	    SP<OS_Interface> interface = new OS_Interface(infile);
	    interface->parser();

	  // initialize the mesh builder and build mesh
	    OS_Builder os_build(interface);
	    mesh = os_build.build_Mesh();

	  // initialize the Opacity builder and build state 
	    Opacity_Builder<OS_Mesh> opacity_build(interface);
	    mat_state = opacity_build.build_Mat(mesh);
	    opacity   = opacity_build.build_Opacity(mesh, mat_state);

	  // do the source initialization
	    sinit = new Source_Init<OS_Mesh>(interface, mesh);
	    sinit->initialize(mesh, opacity, mat_state, rcon, 1);
	}
	
      // sinit diagnostic
	print_sinit(*mesh, *sinit);	
    }
    catch (const assertion &ass)
    {
	cerr << "Dumbass, you screwed up: " << ass.what() << endl;
	return 1;
    }

  // return completed successfully
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tests.cc
//---------------------------------------------------------------------------//
