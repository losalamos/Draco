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
#include "imctest/Source_Init.hh"
#include "imctest/Tally.hh"
#include "imctest/Math.hh"
#include "rng/Rnd_Control.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using IMC::OS_Interface;
using IMC::OS_Builder;
using IMC::OS_Mesh;
using IMC::Mat_State;
using IMC::Opacity_Builder;
using IMC::Opacity;
using IMC::Particle;
using IMC::Tally;
using IMC::Particle_Stack;
using IMC::Source_Init;
using RNG::Sprng;
using RNG::Rnd_Control;
using dsxx::assertion;
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

template<class MT>
void Surface_diagnostic(const MT &mesh)
{
  // input
    string bnd;
    int cell;
    int face;

    cout << "Enter bnd: ";
    cin >> bnd;
    cout << "Enter cell: ";
    cin >> cell;
    cout << "Enter face: ";
    cin >> face;
	
  // check the cells along bnd
    vector<int> surfaces = mesh.get_surcells(bnd);
    cout << surfaces << endl;

  // check the surface face index on bnd in cell
    cout << mesh.get_bndface(bnd, cell) << endl << endl;

  // get the vertices on a cell-face
    typename MT::CCVF_a vertex = mesh.get_vertices(cell, face);
    for (int i = 0; i < vertex[0].size(); i++)
    {
	cout << "(" << vertex[0][i];
	for (int d = 1; d < vertex.size(); d++)
	    cout << "," << vertex[d][i];
	cout << ")" << endl;
    }
}

template<class MT>
void Bank_Particle(const MT &mesh, const Opacity<MT> &xs, Tally<MT> &tally,
		   Rnd_Control &rcon)
{
  // test particle copying and backing

  // make and copy particle

    Sprng ran1 = rcon.get_rn();
    Sprng ran2 = rcon.get_rn();

    vector<double> r1(2), r2(2);
    vector<double> o1(3), o2(3);

    r1[0] = 1.0;
    r1[1] = -1.0;
    r1[2] = 1.0;
    r2 = r1;
    r2[0] = 1.5;
    
    o1[0] = 1.0;
    o2[1] = 1.0;

    Particle<MT> part1(r1, o1, 10.0, mesh.get_cell(r1), ran1);
    Particle<MT> part2(r2, o2, 1.0, mesh.get_cell(r2), ran2);
    Particle<MT> part3(part2);

    Particle<MT, SMrng> part4(part2);

    SP<Particle<MT, SMrng>::Diagnostic> check = 
	new Particle<MT, SMrng>::Diagnostic(cout, true);

    Particle_Stack<Particle<MT> >::Bank sbank;
    sbank.push(part1);
    cout << sbank.size() << endl;
    sbank.push(part2);
    cout << sbank.size() << endl;
    sbank.push(part3);
    cout << sbank.size() << endl;
    sbank.push(part4);
    cout << sbank.size() << endl;

    sbank.pop();
    cout << sbank.size() << endl;
    Particle<MT> part5 = sbank.top();
    sbank.pop();
    cout << sbank.size() << endl;
    cout << part5;	
}

template<class MT>
void write_part(const MT &mesh, Rnd_Control &rcon)
{
  // test io attributes of Particle

  // first make a particle
    vector<double> r(2, 5.0), dir(3, 0.0);
    dir[0] = 1.0;
    Particle<MT> part1(r, dir, 10.0, 1, rcon.get_rn());
    dir[1] = 1.0;
    dir[0] = 0.0;
    r[1] = 3.2;
    Particle<MT> part2(r, dir, 2.1, 32, rcon.get_rn());

    cout << part1 << endl;
    cout << part2 << endl;

  // write particles
    ofstream outfile("part.out");
    part1.write_to_census(outfile);
    part2.write_to_census(outfile);
}

template<class MT>
void read_part(const MT &mesh, Rnd_Control &rcon)
{
  // test io attributes of Particle
    
  // open file and get all the particles
    ifstream infile("part.out", ios::in);

    Particle<MT>::Particle_Buffer buffer(mesh.get_Coord().get_dim());

    while (buffer.read(infile))
    {
	Particle<MT> part(buffer.get_r(), buffer.get_omega(),
			  buffer.get_ew(), buffer.get_cell(), rcon.get_rn(),
			  buffer.get_frac());
	cout << part << endl;
    }

  // delete file
    std::remove("part.out");
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
	    Opacity_Builder<OS_Mesh> opacity_build(interface, mesh);
	    mat_state = opacity_build.build_Mat();
	    opacity   = opacity_build.build_Opacity();

	  // do the source initialization
	    sinit = new Source_Init<OS_Mesh>(interface, mesh);
	    sinit->initialize(mesh, opacity, mat_state, rcon, 1);
	}

      // mesh diagnostics
	Builder_diagnostic(*mesh, *mat_state, *opacity);
	write_part(*mesh, *rcon);
     	read_part(*mesh, *rcon);
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
//                              end of testb.cc
//---------------------------------------------------------------------------//
