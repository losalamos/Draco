//----------------------------------*-C++-*----------------------------------//
// testb.cc
// Thomas M. Evans
// Wed Feb 18 10:45:23 1998
//---------------------------------------------------------------------------//
// @> test driver for Builder classes
//---------------------------------------------------------------------------//

#include "imc/OS_Interface.hh"
#include "imc/OS_Builder.hh"
#include "imc/OS_Mesh.hh"
#include "imc/Mat_State.hh"
#include "imc/Opacity_Builder.hh"
#include "imc/Opacity.hh"
#include "imc/Particle.hh"
#include "imc/Particle_Buffer.hh"
#include "imc/Source_Init.hh"
#include "imc/Tally.hh"
#include "imc/Math.hh"
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
using IMC::Particle_Buffer;
using IMC::Tally;
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
void Bank_Particle(const MT &mesh, const Opacity<MT> &xs, Rnd_Control &rcon)
{
  // test particle copying and backing

  // make and copy particle

    SP<Particle<MT>::Diagnostic> check = 
	new Particle<MT>::Diagnostic(cout, true);

    Particle_Buffer<Particle<MT> >::Bank bank;

    vector<double> r(2, 5.0), dir(3, 0.0);
    Particle<MT> p1(r, dir, 1, 10, rcon.get_rn());
    bank.push(p1);

    cout << endl << "Bank Size should be 1: ";
    cout << bank.size() << endl;
}

template<class MT>
void write_part(const MT &mesh, Rnd_Control &rcon, vector<double> &match)
{
  // test io attributes of Particle
    cout.precision(4);

  // test random number stream
    rcon.set_num(0);
    Sprng ran1 = rcon.get_rn();
    for (int i = 0; i < 25; i++)
	match[i] = ran1.ran();
    
  // first make a particle
    vector<double> r(2, 5.0), dir(3, 0.0);
    dir[0] = 1.0;
    Particle<MT> part1(r, dir, 10.0, 1, ran1);
    dir[1] = 1.0;
    dir[0] = 0.0;
    r[1] = 3.2;
    Particle<MT> part2(r, dir, 2.1, 32, rcon.get_rn());

    cout << part1 << endl;
    cout << part2 << endl;

  // write particles
    ofstream outfile("part.out");
    Particle_Buffer<Particle<MT> > buffer(mesh, rcon);
    buffer.write_census(outfile, part1);
    buffer.write_census(outfile, part2);
}

template<class MT>
void read_part(const MT &mesh, const Rnd_Control &rcon,
	       vector<double> &match)
{
  // test io attributes of Particle
    
    cout.precision(4);

  // open file and get all the particles
    ifstream infile("part.out", ios::in);

  // read particles

    Particle_Buffer<Particle<MT> > buffer(mesh, rcon);
    SP<Particle_Buffer<Particle<MT> >::Census_Buffer> cenpart;
    int index = 0;
    do
    {
	index++;
	cenpart = buffer.read_census(infile);
	if (cenpart)
	{
	    if (index == 1)
		for (int i = 25; i < 50; i++)
		    match[i] = cenpart->random.ran();
	    Particle<MT> part(cenpart->r, cenpart->omega, cenpart->ew,
			      cenpart->cell, cenpart->random,
			      cenpart->fraction);
	    cout << part << endl;
	}
    } while (cenpart);
	
  // delete file
    std::remove("part.out");
}

void print_ran(vector<double> &control, Rnd_Control &rcon)
{
    cout.precision(4);

    rcon.set_num(0);
    Sprng random = rcon.get_rn();

    for (int i = 0; i < 50; i++)
	control[i] = random.ran();
}

template<class MT>
void Persistence_diagnostic(const MT &mesh, Rnd_Control &rcont)
{
    vector<double> control(50);
    vector<double> match(50);

    write_part(mesh, rcont, match);
    read_part(mesh, rcont, match);
    print_ran(control, rcont);

    cout << "The following should match:" << endl;
    cout.precision(4);
    for (int i = 0; i < 50; i++)
	cout << setw(5) << i+1 << setw(10) << control[i] 
	     << setw(10) << match[i] << endl;
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
	Persistence_diagnostic(*mesh, *rcon);
	Bank_Particle(*mesh, *opacity, *rcon);
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
