//----------------------------------*-C++-*----------------------------------//
// testp.cc
// Thomas M. Evans
// Mon Apr 13 17:31:21 1998
//---------------------------------------------------------------------------//
// @> test executable to try out parallelism in IMCTEST
//---------------------------------------------------------------------------//

#include "imctest/OS_Interface.hh"
#include "imctest/OS_Builder.hh"
#include "imctest/OS_Mesh.hh"
#include "imctest/Mat_State.hh"
#include "imctest/Opacity_Builder.hh"
#include "imctest/Opacity.hh"
#include "imctest/Parallel_Builder.hh"
#include "ds++/SP.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>

using IMC::OS_Interface;
using IMC::OS_Builder;
using IMC::OS_Mesh;
using IMC::Mat_State;
using IMC::Opacity_Builder;
using IMC::Opacity;
using IMC::Parallel_Builder;
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
    cout << "Wrote Mesh diagnostic file " << title << " from proc. " 
	 << mynode << endl;
}

int main(int argc, char *argv[])
{

  // init C4 stuff
    Init(argc, argv);
    mynode  = C4::node();
    mynodes = C4::nodes();

  // declare geometry and material stuff
    SP<OS_Mesh> mesh;
    SP< Mat_State<OS_Mesh> > mat_state;
    SP< Opacity<OS_Mesh> > opacity;

  // read input and stuff on the host-topology
    if (!mynode)
    {
	string infile = argv[1];
	
      // run the interface parser
	SP<OS_Interface> interface = new OS_Interface(infile);
	interface->parser();

      // initialize the mesh builder and build mesh
	OS_Builder os_build(interface);
	mesh = os_build.build_Mesh();

      // initialize the Opacity builder and build state 
	Opacity_Builder<OS_Mesh> opacity_build(interface, mesh);
	mat_state = opacity_build.build_Mat();
	opacity   = opacity_build.build_Opacity();
    }

  // make parallel builder object to do my mesh decomposition
    Parallel_Builder<OS_Mesh> pcomm;
    
    if (!mynode)
	pcomm.send_Mesh(*mesh);

    if (mynode)
    {
        mesh = pcomm.recv_Mesh();
	std::cout << mesh->num_cells() << std::endl;
    } 
    
    if (mesh) 
	Builder_diagnostic(*mesh, *mat_state, *opacity);

  // c4 end
    Finalize();
}

//---------------------------------------------------------------------------//
//                              end of testp.cc
//---------------------------------------------------------------------------//
