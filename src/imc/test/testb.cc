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
#include "imctest/Random.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

using IMC::OS_Interface;
using IMC::OS_Builder;
using IMC::OS_Mesh;
using IMC::Mat_State;
using IMC::Opacity_Builder;
using IMC::Opacity;
using IMC::Particle;
using IMC::Random;
using namespace std;

template<class MT>
void Builder_diagnostic(const MT &mesh, const Mat_State<MT> &mat,
			const Opacity<MT> &opacity)
{
  // do some diagnostic checks

    ofstream output("try.dat");

    output << "Coordinate System: " << mesh.Coord().Get_coord() << endl;
    output << "Mesh Size: " << mesh.Num_cells() << endl;
    output << endl;

    for (int cell = 1; cell <= mesh.Num_cells(); cell++)
	mesh.Print(output, cell);
    output << endl;

    output << mat;
    output << endl;

    output << opacity;
    output << endl;

    output << "Mesh:      " << mesh.Num_cells() << endl;
    output << "Opacity:   " << opacity.Num_cells() << endl;
    output << "Mat_State: " << mat.Num_cells() << endl;
}

template<class MT>
void Run_Particle(const MT &mesh, const Opacity<MT> &opacity, long seed)
{
  // transport a particle

  // set diagnostics
    
    ofstream output("history", ios::app);
    SP<Particle<MT>::Diagnostic> check = 
	new Particle<MT>::Diagnostic(output, true);

  // initialize particle
    Particle<MT> particle(mesh, seed, 1.0, 10.0);
    assert (particle.Status());

  // origin and source
    
    vector<double> origin(mesh.Coord().Get_dim());   
    vector<double> source(mesh.Coord().Get_sdim());

    for (int i = 0; i < origin.size(); i++)
    {
	cout << "Origin " << i+1 << ": ";
	cin >> origin[i];
    }
    for (int i = 0; i < source.size(); i++)
    {
	cout << "Omega " << i+1 << ": ";
	cin >> source[i];
    }
    cout << endl;
    
  // source
    particle.Source(origin, source, mesh);

  // transport
    particle.Transport(mesh, opacity, check);

  // assert that particle is indeed dead
    assert (!particle.Status());
}

main()
{
  // declare geometry and material stuff
    SP<OS_Mesh> mesh;
    SP< Mat_State<OS_Mesh> > mat_state;
    SP< Opacity<OS_Mesh> > opacity;

  // scoping blocks for build-stuff
    {
	string infile;
	cout << "Name the input file" << endl;
	cin >> infile;

      // run the interface parser
	SP<OS_Interface> interface = new OS_Interface(infile);
	interface->Parser();

      // initialize the mesh builder and build mesh
	OS_Builder os_build(interface);
	mesh = os_build.Build_Mesh();

      // initialize the Opacity builder and build state 
	Opacity_Builder<OS_Mesh> opacity_build(interface, mesh);
	mat_state = opacity_build.Build_Mat();
	opacity   = opacity_build.Build_Opacity();
    }

  // mesh diagnostics
    Builder_diagnostic(*mesh, *mat_state, *opacity);

    for (int cell = 1; cell <= mesh->Num_cells(); cell++)
	cout << cell << " " << mesh->Volume(cell) << endl;

//     long seed = -345632;
//     Run_Particle(*mesh, *opacity, seed);
}


//---------------------------------------------------------------------------//
//                              end of testb.cc
//---------------------------------------------------------------------------//
