//----------------------------------*-C++-*----------------------------------//
// testa.cc
// Thomas M. Evans
// Tue Jun  9 18:42:16 1998
//---------------------------------------------------------------------------//
// @> test parallel building and mixing for IMC
//---------------------------------------------------------------------------//

#include "imc/OS_Interface.hh"
#include "imc/OS_Builder.hh"
#include "imc/OS_Mesh.hh"
#include "imc/Mat_State.hh"
#include "imc/Opacity_Builder.hh"
#include "imc/Opacity.hh"
#include "imc/Parallel_Builder.hh"
#include "imc/Source_Init.hh"
#include "imc/Particle_Buffer.hh"
#include "imc/Particle.hh"
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

using IMC::OS_Interface;
using IMC::OS_Builder;
using IMC::OS_Mesh;
using IMC::Mat_State;
using IMC::Opacity_Builder;
using IMC::Opacity;
using IMC::Parallel_Builder;
using IMC::Source_Init;
using IMC::Particle_Buffer;
using IMC::Particle;
using IMC::Global::rn_stream;
using RNG::Rnd_Control;
using RNG::Sprng;
using namespace std;
using namespace C4;

// declare node
int mynode;
int mynodes;

template<class MT>
void topology(const MT &mesh, const Parallel_Builder<MT> &pcom)
{
    cout << ">> IMC Topology Diagnostic" << endl << endl;

    cout << "Number of procs: " << nodes << endl;

    for (int i = 0; i < nodes(); i++)
    {
	cout << setw(9) << pcom.num_cells(i) 
	     << " Cells on node " << setw(5) << i << endl;
	cout << "-----------------------------" << endl;
	cout << "          local        global" << endl;
	for (int j = 1; j <= pcom.num_cells(i); j++)
	{
	    cout << setw(15) << j << setw(14) << pcom.master_cell(j,i)
		 << endl;
	}
	cout << endl;
    }

    for (int i = 1; i <= mesh.num_cells(); i++)
    {
	cout << "Global Cell " << i << " on " << pcom.num_procs(i)
	     << " procs" << endl;
	cout << "      Proc" << "     Local" << endl;
	for (int p = 0; p < nodes(); p++)
	    if (pcom.imc_cell(i, p))
		cout << setw(10) << p << setw(10) 
		     << pcom.imc_cell(i, p) << endl;
    }
    cout << endl;
}

int main(int argc, char *argv[])
{    
  // init C4 stuff
    C4::Init(argc, argv);
    mynode  = C4::node();
    mynodes = C4::nodes();
 
 // try block
    try
    {
      // lets look at our buffers
	int sb = 1000;
	int db = 1000 * (9);
	int ib = 1000 * (2);
	int cb = 1000 * (500);
	Check (IMC::Global::buffer_s == sb);
	Check (IMC::Global::buffer_d == db);
	Check (IMC::Global::buffer_i == ib);
	Check (IMC::Global::buffer_c == cb);

      // declare geometry and material stuff
	SP<OS_Mesh> mesh;
	SP< Mat_State<OS_Mesh> > mat_state;
	SP< Opacity<OS_Mesh> > opacity;
	SP< Source_Init<OS_Mesh> > sinit;
	SP<Rnd_Control> rcon = new Rnd_Control(9836592);

      // read input and stuff on the host-topology
	if (!mynode)
	{
	    cout << ">> Running through problem builders" << endl;
	    string infile = argv[1];
	
	  // run the interface parser
	    SP<OS_Interface> interface = new OS_Interface(infile);
	    interface->parser();
            cout << "** Read input file on host " << mynode << endl;

	  // initialize the mesh builder and build mesh
	    OS_Builder os_build(interface);
	    mesh = os_build.build_Mesh();
            cout << "** Built mesh on host " << mynode << endl;

	  // initialize the Opacity builder and build state 
	    Opacity_Builder<OS_Mesh> opacity_build(interface, mesh);
	    mat_state = opacity_build.build_Mat();
	    opacity   = opacity_build.build_Opacity();
            cout << "** Built opacities on host " << mynode << endl;

	  // do the source initialization
	    sinit = new Source_Init<OS_Mesh>(interface, mesh);
	    sinit->initialize(mesh, opacity, mat_state, rcon, 1);
            cout << "** Initialized source on host " << mynode << endl;
            cout << endl;
	}

      // make parallel builder object to do send/receives of objects
	SP< Parallel_Builder<OS_Mesh> > pcomm;
	SP< Particle_Buffer<Particle<OS_Mesh> > > buffer;
    
	if (!mynode)
	{
	  // make parallel builder object to do my mesh decomposition
 	    pcomm = new Parallel_Builder<OS_Mesh>(*mesh, *sinit);

	  // topology diagnostic
	    topology(*mesh, *pcomm);
	
          // make Particle buffer
            buffer = new Particle_Buffer<Particle<OS_Mesh> >(*mesh, *rcon);

          // send out objects
            pcomm->send_Mesh(*mesh);
	    pcomm->send_Opacity(*opacity);
	}
	
  	if (mynode)
  	{	
  	  // make parallel builder object to receive objects
  	    pcomm = new Parallel_Builder<OS_Mesh>();

          // get mesh and opacity
  	    mesh    = pcomm->recv_Mesh();
  	    opacity = pcomm->recv_Opacity(mesh);

          // make a particle buffer on this node
            buffer  = new Particle_Buffer<Particle<OS_Mesh> >(*mesh, *rcon);

          // get the source for this node
  	}

      // lets look at our buffers
	Check (IMC::Global::buffer_s == sb);
	Check (IMC::Global::buffer_d == db);
	Check (IMC::Global::buffer_i == ib);
	Check (IMC::Global::buffer_c == cb);
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
