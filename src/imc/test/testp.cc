//----------------------------------*-C++-*----------------------------------//
// testp.cc
// Thomas M. Evans
// Mon Apr 13 17:31:21 1998
//---------------------------------------------------------------------------//
// @> test executable to try out parallelism in IMCTEST
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
    cout << "** Wrote Full Mesh diagnostic file " << title << " from proc. " 
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
    cout << "** Wrote Mesh diagnostic file " << title << " from proc. " 
	 << mynode << endl;
}

template<class MT, class PT>
void comm_particle1(const MT &mesh,
		    const Particle_Buffer<PT> &buffer, 
		    Rnd_Control &rcon)
{
  // lets do stuff on the host
    if (!mynode)
    {    
	cout << endl << ">> Particle Comm_Buffer blocking send/recv test 1"
	     << endl;

      // let's make some particles
	vector<double> r(mesh.get_Coord().get_dim(), 1.0);
	vector<double> omega(mesh.get_Coord().get_sdim(), 0.0);
	Sprng ran1 = rcon.get_rn();
	Sprng ran2 = rcon.get_rn();
	Sprng ran3 = rcon.get_rn();
	for (int i = 1; i <= 10; i++) ran1.ran();
	for (int i = 1; i <= 10; i++) ran2.ran();
	for (int i = 1; i <= 10; i++) ran3.ran();
	PT p1(r, omega, 1.0, 1, ran1);
	fill(r.begin(), r.end(), 2.0);
	PT p2(r, omega, 2.0, 2, ran2);
	fill(r.begin(), r.end(), 3.0);
	PT p3(r, omega, 3.0, 3, ran3);

      // lets dump them
	cout << "** Printing Particles on node " << mynode << endl;
	cout << p1 << p2 << p3 << endl;

      // lets buffer them
	cout << "** Putting Particles in a Comm_Buffer on node " 
	     << mynode << endl;
	Particle_Buffer<PT>::Comm_Buffer cbuf;
	buffer.buffer_particle(cbuf, p1);
	buffer.buffer_particle(cbuf, p2);
 	buffer.buffer_particle(cbuf, p3);

      // lets send them
	for (int np = 1; np < mynodes; np++) 
	{
	    buffer.send_buffer(cbuf, np);
	    cout << "** Sent buffer to node " << np << endl;
	}
    }

  // lets get stuff from the host
    if (mynode)
    {
	SP<Particle_Buffer<PT>::Comm_Buffer> rbuf;
	
      // let's receive the buffers
	rbuf = buffer.recv_buffer(0);
	cout << "** We received a buffer on " << mynode << " with " 
	     << rbuf->n_part << " particles from node " << 0 << endl;

      // write the buffers to census files
	ostringstream stitle;
	stitle << "census." << mynode;
	string title = stitle.str();
	
	ofstream cenout(title.c_str());
	buffer.write_census(cenout, *rbuf);
	cout << "** Wrote the file " << title << " on node " << mynode
	     << endl;
    }

  // lets open the census files and see what we've got
    HTSyncSpinLock h;
    {
	ostringstream stitle;
	stitle << "census." << mynode;
	string title = stitle.str();

	ifstream infile(title.c_str());

	if (!infile == false)
	{
	    cout << "** Opening census file " << title << " on node " 
		 << mynode << endl;

	  // read particles
	    SP<Particle_Buffer<PT>::Census_Buffer> cenpart;
	    cout << "** Printing Particles on node " << mynode 
		 << endl;
	    do
	    {
		cenpart = buffer.read_census(infile);
		if (cenpart)
		{
		    PT part(cenpart->r, cenpart->omega, cenpart->ew,
			    cenpart->cell, cenpart->random,
			    cenpart->fraction);
		    cout << part << endl << endl;

		  // random number check
		    rcon.set_num(cenpart->random.get_num());
		    Sprng rand = rcon.get_rn();
		    cout << "The following should match:" << endl;
		    cout << "---------------------------" << endl;
		    for (int i = 1; i <= 20; i++)
		    {
			if (i <= 10)
			    rand.ran();
			else
			{
			    cout.precision(4);
			    cout << setw(10) << rand.ran() << setw(10)
				 << cenpart->random.ran() << endl;
			}
		    }
		    cout << "---------------------------" << endl << endl;
		}
	    } while (cenpart);

	  // closing Census file
	    remove(title.c_str());
	    cout << "** Removing file " << title << " on node " << mynode 
		 << endl;
	}  
    }	    
}

template<class MT, class PT>
void comm_particle2(const MT &mesh,
		    const Particle_Buffer<PT> &buffer, 
		    Rnd_Control &rcon)
{
  // lets do stuff on the host
    if (!mynode)
    {    
	cout << endl << ">> Particle Comm_Buffer blocking send/recv test 2"
	     << endl;

      // let's make some particles
	vector<double> r(mesh.get_Coord().get_dim(), 10.0);
	vector<double> omega(mesh.get_Coord().get_sdim(), 0.0);
	Sprng ran1 = rcon.get_rn();
	Sprng ran2 = rcon.get_rn();
	Sprng ran3 = rcon.get_rn();
	for (int i = 1; i <= 10; i++) ran1.ran();
	for (int i = 1; i <= 10; i++) ran2.ran();
	for (int i = 1; i <= 10; i++) ran3.ran();
	PT p1(r, omega, 10.0, 10, ran1);
	fill(r.begin(), r.end(), 20.0);
	PT p2(r, omega, 20.0, 20, ran2);
	fill(r.begin(), r.end(), 30.0);
	PT p3(r, omega, 30.0, 30, ran3);

      // lets dump them
	cout << "** Printing Particles on node " << mynode << endl;
	cout << p1 << p2 << p3 << endl;

      // lets buffer them
	cout << "** Putting Particles in a Comm_Buffer on node " 
	     << mynode << endl;
	Particle_Buffer<PT>::Comm_Buffer cbuf;
	buffer.buffer_particle(cbuf, p1);
	buffer.buffer_particle(cbuf, p2);
 	buffer.buffer_particle(cbuf, p3);

      // lets send them
	for (int np = 1; np < mynodes; np++) 
	{
	    buffer.send_buffer(cbuf, np);
	    cout << "** Sent buffer to node " << np << endl;
	}
    }

  // lets get stuff from the host
    if (mynode)
    {
	SP<Particle_Buffer<PT>::Comm_Buffer> rbuf;
	
      // let's receive the buffers
	rbuf = buffer.recv_buffer(0);
	cout << "** We received a buffer on " << mynode << " with " 
	     << rbuf->n_part << " particles from node " << 0 << endl;

      // let's put them into a bank
	Particle_Buffer<PT>::Comm_Bank bank;
	buffer.add_to_bank(*rbuf, bank);
	cout << "** We made a Particle Bank on node " << mynode << endl;

      // let's get the guys out of the bank and check'em out
	while (bank.size() > 0)
	{
	  // print out the particle
	    PT &particle = bank.top();
	    cout << "** Getting particle out of bank on node " << mynode
		 << endl;
	    cout << particle << endl;

	  // do a random number check
	    rcon.set_num(particle.get_random().get_num());
	    Sprng rand = rcon.get_rn();
	    cout << "The following should match:" << endl;
	    cout << "---------------------------" << endl;
	    for (int i = 1; i <= 20; i++)
	    {
		if (i <= 10)
		    rand.ran();
		else
		{
		    cout.precision(4);
		    cout << setw(10) << rand.ran() << setw(10)
			 << particle.get_random().ran() << endl;
		}
	    }
	    cout << "---------------------------" << endl << endl;
	    
	    bank.pop();
	    cout << "** There are now " << bank.size() 
		 << " particles left in the bank on node " << mynode
		 << endl;
	}
    }
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
	Check (Particle_Buffer<Particle<OS_Mesh> >::get_buffer_s() == sb);
	Check (Particle_Buffer<Particle<OS_Mesh> >::get_buffer_d() == db);
	Check (Particle_Buffer<Particle<OS_Mesh> >::get_buffer_i() == ib);
	Check (Particle_Buffer<Particle<OS_Mesh> >::get_buffer_c() == cb);

	Particle_Buffer<Particle<OS_Mesh> >::set_buffer_size(3);

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

          // make Particle buffer
            buffer  = new Particle_Buffer<Particle<OS_Mesh> >(*mesh, *rcon);

          // send out objects
            pcomm->send_Mesh(*mesh);
	    pcomm->send_Opacity(*opacity);
	    Builder_diagnostic(*mesh, *mat_state, *opacity);
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
     
  	if (mesh) 
	{
  	    Builder_diagnostic(*mesh, *opacity);
	    comm_particle1(*mesh, *buffer, *rcon);
	    comm_particle2(*mesh, *buffer, *rcon);
	}

      // print out the buffers
	HTSyncSpinLock h;
	{
	    cout << endl << ">> Buffer sizes" << endl;
	    cout << "** On node " << mynode << " the buffer size was "
		 << Particle_Buffer<Particle<OS_Mesh> >::get_buffer_s() 
		 << endl;
	    cout << "** Buffer_d = " 
		 << Particle_Buffer<Particle<OS_Mesh> >::get_buffer_d() 
		 << endl;
	    cout << "** Buffer_i = " 
		 << Particle_Buffer<Particle<OS_Mesh> >::get_buffer_i() 
		 << endl;
	    cout << "** Buffer_c = " 
		 << Particle_Buffer<Particle<OS_Mesh> >::get_buffer_c() 
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
//                              end of testp.cc
//---------------------------------------------------------------------------//
