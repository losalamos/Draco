//----------------------------------*-C++-*----------------------------------//
// testb.cc
// Thomas M. Evans
// Wed Feb 18 10:45:23 1998
//---------------------------------------------------------------------------//
// @> test driver for Builder classes
//---------------------------------------------------------------------------//

#include "imctest/OS_Interface.hh"
#include "imctest/OS_Builder.hh"
#include "imctest/OS_Mesh.hh"
#include "imctest/Mat_State.hh"
#include "imctest/Opacity_Builder.hh"
#include "imctest/Opacity.hh"
#include "imctest/Particle.hh"
#include "imctest/Tally.hh"
#include "imctest/Random.hh"
#include "imctest/Math.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <fstream>
#include <string>

using IMC::OS_Interface;
using IMC::OS_Builder;
using IMC::OS_Mesh;
using IMC::Mat_State;
using IMC::Opacity_Builder;
using IMC::Opacity;
using IMC::Particle;
using IMC::Tally;
using IMC::Random;
using IMC::Particle_Stack;
using IMC::Global::operator<<;
using namespace std;

template<class MT>
void Builder_diagnostic(const MT &mesh, const Mat_State<MT> &mat,
			const Opacity<MT> &opacity)
{
  // do some diagnostic checks

    string title = "try.dat";

    ofstream output(title.c_str());

    output << "Coordinate System: " << mesh.get_Coord().get_Coord() << endl;
    output << "Mesh Size: " << mesh.num_cells() << endl;
    output << endl;

    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	mesh.print(output, cell);
    output << endl;

    output << mat;
    output << endl;

    output << opacity;
    output << endl;

    output << "Mesh:      " << mesh.num_cells() << endl;
    output << "Opacity:   " << opacity.num_cells() << endl;
    output << "Mat_State: " << mat.num_cells() << endl;
}

int main(int argc, char *argv[])
{
  // declare geometry and material stuff
    SP<OS_Mesh> mesh;
    SP< Mat_State<OS_Mesh> > mat_state;
    SP< Opacity<OS_Mesh> > opacity;
    
  // scoping blocks for build-stuff
    {
	string infile = "in2";
	
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

  // do a diagnostic on Mesh build

    Builder_diagnostic(*mesh, *mat_state, *opacity);
}


//---------------------------------------------------------------------------//
//                              end of testb.cc
//---------------------------------------------------------------------------//
