//----------------------------------*-C++-*----------------------------------//
// testp.cc
// Thomas M. Evans
// Mon Apr 13 17:31:21 1998
//---------------------------------------------------------------------------//
// @> test executable to try out parallelism in IMCTEST.
//---------------------------------------------------------------------------//

#include "imctest/OS_Interface.hh"
#include "imctest/OS_Builder.hh"
#include "imctest/OS_Mesh.hh"
#include "imctest/Mat_State.hh"
#include "imctest/Opacity_Builder.hh"
#include "imctest/Opacity.hh"
#include "imctest/Parallel_Builder.hh"
#include "imctest/Source_Init.hh"
#include "rng/Rnd_Control.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using IMC::OS_Interface;
using IMC::OS_Builder;
using IMC::OS_Mesh;
using IMC::Mat_State;
using IMC::Opacity_Builder;
using IMC::Opacity;
using IMC::Parallel_Builder;
using IMC::Source_Init;
using RNG::Rnd_Control;
using namespace std;
using namespace C4;

// declare node
int mynode;
int mynodes;

template<class MT>
void Builder_diagnostic(const MT &mesh, const Mat_State<MT> &mat,
			const Opacity<MT> &opacity)
{
  // do some diagnostic checks

    string title = "full.dat";
    ofstream output(title.c_str());

  // title header
    output << "Coordinate System: " << mesh.get_Coord().get_Coord() << endl;
    output << "Mesh Size: " << mesh.num_cells() << endl;
    output << endl;

  // print mesh
    output << mesh;
    output << endl;

  // print mat state
    output << mat;
    output << endl;

  // print opacity
    output << opacity;
    output << endl;

  // final diagnostics
    output << "Mesh:      " << mesh.num_cells() << endl;
    output << "Opacity:   " << opacity.num_cells() << endl;
    output << "Mat_State: " << mat.num_cells() << endl;

  // message
    cout << "Wrote Full Mesh diagnostic file " << title << " from proc. " 
	 << mynode << endl;
}

template<class MT>
void Builder_diagnostic(const MT &mesh,	const Opacity<MT> &opacity)
{
  // do some diagnostic checks

    ostringstream stitle;
    stitle << mynode << ".dat";
    string title = stitle.str();

    ofstream output(title.c_str());

  // title header
    output << "Coordinate System: " << mesh.get_Coord().get_Coord() << endl;
    output << "Mesh Size: " << mesh.num_cells() << endl;
    output << endl;

  // print mesh
    output << mesh;
    output << endl;

  // print opacity
    output << opacity;
    output << endl;

  // final diagnostics
    output << "Mesh:      " << mesh.num_cells() << endl;
    output << "Opacity:   " << opacity.num_cells() << endl;

  // message
    cout << "Wrote Mesh diagnostic file " << title << " from proc. " 
	 << mynode << endl;
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
      // declare geometry and material stuff
	SP<OS_Mesh> mesh;
	SP< Mat_State<OS_Mesh> > mat_state;
	SP< Opacity<OS_Mesh> > opacity;
	SP< Source_Init<OS_Mesh> > sinit;
	SP<Rnd_Control> rcon = new Rnd_Control(9836592);

      // read input and stuff on the host-topology
	if (!mynode)
	{
	    string infile = argv[1];
	
	  // run the interface parser
	    SP<OS_Interface> interface = new OS_Interface(infile);
	    interface->parser();
            cout << "Read input file on host " << mynode << endl;

	  // initialize the mesh builder and build mesh
	    OS_Builder os_build(interface);
	    mesh = os_build.build_Mesh();
            cout << "Built mesh on host " << mynode << endl;

	  // initialize the Opacity builder and build state 
	    Opacity_Builder<OS_Mesh> opacity_build(interface, mesh);
	    mat_state = opacity_build.build_Mat();
	    opacity   = opacity_build.build_Opacity();
            cout << "Built opacities on host " << mynode << endl;

	  // do the source initialization
	    sinit = new Source_Init<OS_Mesh>(interface, mesh);
	    sinit->initialize(mesh, opacity, mat_state, rcon, 1);
            cout << "Initialized source on host " << mynode << endl;
            cout << endl;
	}

      // make parallel builder object to do send/receives of objects
	SP<Parallel_Builder<OS_Mesh> > pcomm;
    
	if (!mynode)
	{
	  // make parallel builder object to do my mesh decomposition
	  // 	    pcomm = new Parallel_Builder<OS_Mesh>(*mesh, *sinit);
 	    pcomm = new Parallel_Builder<OS_Mesh>();
	    pcomm->parallel_params(*mesh, *sinit);
	    Check (pcomm);
	    Check (mesh);
	    Check (opacity);
            pcomm->send_Mesh(*mesh);
	    pcomm->send_Opacity(*opacity);
	    Builder_diagnostic(*mesh, *mat_state, *opacity);
	    
	  // print out salient quantities from Parallel_Builder
//          for (int i = 0; i < nodes(); i++)
//              cout << setw(10) << i << setw(10) 
//                   << pcomm->cells_per_proc[i].size() << endl;
	}
	
  	if (mynode)
  	{	
  	  // make parallel builder object to receive objects
  	    pcomm = new Parallel_Builder<OS_Mesh>();
	    Check (pcomm);
  	    mesh    = pcomm->recv_Mesh();
  	    opacity = pcomm->recv_Opacity(mesh);
	    Check (mesh);
	    Check (opacity);
  	}
     
  	if (mesh) 
  	    Builder_diagnostic(*mesh, *opacity);

    }
    catch (const dsxx::assertion &ass)
    {
	cout << "Dumbass, you screwed up: " << ass.what() << endl;
	return 1;
    }
//     catch(...)
//     {
// 	cout << "HELP ME" << endl;
// 	return 1;
//     }

  // c4 end
    C4::Finalize();

  // we ran ok
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of testp.cc
//---------------------------------------------------------------------------//
