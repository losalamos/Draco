//----------------------------------*-C++-*----------------------------------//
// testp.cc
// Thomas M. Evans
// Mon Apr 13 17:31:21 1998
//---------------------------------------------------------------------------//
// @> test executable to try out parallelism in IMCTEST
//---------------------------------------------------------------------------//

#include "imc/OS_Mesh.hh"
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

using IMC::OS_Mesh;
using IMC::Particle_Buffer;
using IMC::Particle;
using RNG::Rnd_Control;
using RNG::Sprng;
using dsxx::SP;
using namespace std;
using namespace C4;

// declare node
int mynode;
int mynodes;

void comp_comm()
{
    Particle_Buffer<Particle<OS_Mesh> >::set_buffer_size(3);

  // declare communication
    SP<OS_Mesh> mesh; 
    SP<Rnd_Control> rcon = new Rnd_Control(9836592);
    Particle_Buffer<Particle<OS_Mesh> > buffer(8, 2, rcon->get_size()); 
    Particle_Buffer<Particle<OS_Mesh> >::Comm_Buffer send_buffer;
    Particle_Buffer<Particle<OS_Mesh> >::Comm_Buffer recv_buffer;
    Particle_Buffer<Particle<OS_Mesh> >::Bank bank;

  // post a receive on each processor
    if (node() == 0)
	buffer.post_arecv(recv_buffer, 1);
    if (node() == 1)
	buffer.post_arecv(recv_buffer, 0);

  // make some particles
    {
	vector<double> r(3, static_cast<double>(node()));
	vector<double> o(3, 5.0);
	Particle<OS_Mesh>  part(r, o, 10.0 * (node() + 1), 50, rcon->get_rn());
	buffer.buffer_particle(send_buffer, part);
    }

  // async send the particles
    if (node() == 0)
    {
	for (int i = 0; i < 100000; i++);
	buffer.asend_buffer(send_buffer, 1);
    }
    if (node() == 1)
	buffer.asend_buffer(send_buffer, 0);

  // get the particles 
    bool done = false;
    while (!done)
	done = buffer.async_check(recv_buffer);
    buffer.add_to_bank(recv_buffer, bank);

  // check to see we got the right stuff
    while (bank.size())
    {
	if (!node())
	{
	    vector<double> r(3, 1.0);
	    vector<double> o(3, 5.0);
	    Particle<OS_Mesh> part(r, o, 20.0, 50, rcon->get_rn(0));
	    if (*bank.top() == part)
		cout << "TEST IS SUCCESSFUL ON NODE " << node() << endl;
            else
                cout << "TEST IS UNSUCCESSFUL ON NODE " << node() << endl;	    
	    bank.pop();
	}
	if (node())
	{
	    vector<double> r(3, 0.0);
	    vector<double> o(3, 5.0);
	    Particle<OS_Mesh> part(r, o, 10.0, 50, rcon->get_rn(0));
	    if (*bank.top() == part)
		cout << "TEST IS SUCCESSFUL ON NODE " << node() << endl;
            else
                cout << "TEST IS UNSUCCESSFUL ON NODE " << node() << endl;	    
	    bank.pop();
	}
    }

  // free the recv buffers
    buffer.async_free(recv_buffer);
}


int main(int argc, char *argv[])
{    
  // init C4 stuff
    C4::Init(argc, argv);
    mynode  = C4::node();
    mynodes = C4::nodes();

  // tests
    comp_comm();

  // c4 end
    C4::Finalize();

  // we ran ok
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of testp.cc
//---------------------------------------------------------------------------//
