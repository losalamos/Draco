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
#include "imctest/Math.hh"
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
using IMC::Global::operator<<;
using namespace std;

template<class MT>
void Builder_diagnostic(const MT &mesh, const Mat_State<MT> &mat,
			const Opacity<MT> &opacity)
{
  // do some diagnostic checks

    ofstream output("try.dat");

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

template<class MT>
void Run_Particle(const MT &mesh, const Opacity<MT> &opacity, long seed)
{
  // transport a particle

  // set diagnostics
    
    ofstream output("history", ios::app);
    SP<Particle<MT>::Diagnostic> check = 
	new Particle<MT>::Diagnostic(output, true);

  // initialize particle
    Particle<MT> particle(mesh, seed, 1.0);
    assert (particle.status());

  // origin and source
    
    vector<double> origin(mesh.get_Coord().get_dim());   
    vector<double> source(mesh.get_Coord().get_sdim());

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
    particle.source(origin, source, mesh);

  // transport
    particle.transport_IMC(mesh, opacity, check);

  // assert that particle is indeed dead
    assert (!particle.status());
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
	interface->parser();

      // initialize the mesh builder and build mesh
	OS_Builder os_build(interface);
	mesh = os_build.build_Mesh();

      // initialize the Opacity builder and build state 
	Opacity_Builder<OS_Mesh> opacity_build(interface, mesh);
	mat_state = opacity_build.build_Mat();
	opacity   = opacity_build.build_Opacity();
    }

  // mesh diagnostics
    Builder_diagnostic(*mesh, *mat_state, *opacity);

    vector<int> surfaces = mesh->get_surcells("lox");
    cout << surfaces << endl;
    surfaces = mesh->get_surcells("hix");
    cout << surfaces << endl;  
    surfaces = mesh->get_surcells("loy");
    cout << surfaces << endl;
    surfaces = mesh->get_surcells("hiy");
    cout << surfaces << endl;
    surfaces = mesh->get_surcells("loz");
    cout << surfaces << endl;  
    surfaces = mesh->get_surcells("hiz");
    cout << surfaces << endl;

    cout << mesh->get_bndface("hiz", 9);

//     long seed = -345632;
//     Run_Particle(*mesh, *opacity, seed);
}


//---------------------------------------------------------------------------//
//                              end of testb.cc
//---------------------------------------------------------------------------//
