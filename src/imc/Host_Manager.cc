//----------------------------------*-C++-*----------------------------------//
// Host_Manager.cc
// Thomas M. Evans
// Sun Aug  2 12:43:05 1998
//---------------------------------------------------------------------------//
// @> Host_Manager class implementation file
//---------------------------------------------------------------------------//

#include "imc/Host_Manager.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

IMCSPACE

// draco components
using C4::node;
using C4::nodes;
using C4::Send;
using C4::Recv;
using C4::C4_Req;
using C4::gsync;
using C4::RecvAsync;
using C4::SendAsync;

// stl components
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::ios;
using std::setiosflags;
using std::setw;
using std::fill;
using std::fabs;
using std::vector;
using std::string;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor

template<class MT, class BT, class IT, class PT>
Host_Manager<MT,BT,IT,PT>::Host_Manager(int cycle_)
    : delta_t(0), cycle(cycle_), dump_f(Global::huge)
{
  // all the SPs should be null defined

  // objects used by all processors
    Check (!mesh);
    Check (!opacity);
    Check (!mat_state);
    Check (!rnd_con);
    Check (!buffer);
    Check (!source);
    Check (!tally);
    Check (!communicator);

  // initialization objects
    Check (!source_init);

  // objects used to control the census
    Check (!new_census_bank);
}

//---------------------------------------------------------------------------//
// execute IMC in the host
//---------------------------------------------------------------------------//
// run an IMC cycle

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::execute_IMC(typename IT::Arguments &arg)
{
  // run through an IMC timestep
  
  // initialize IMC
    initialize(arg);
    
  // return if we are on cycle 0
    if (cycle == 0)
	return;
    
  // do a time-step
    step_IMC();

  // cleanup before shipping off
    regroup(arg);
}

//---------------------------------------------------------------------------//
// IMC initialization
//---------------------------------------------------------------------------//
// here we build the Mesh, Opacity, Mat_State, etc. ie. everything needed on
// each processor to do IMC transport

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::initialize(const typename IT::Arguments &arg) 
{
  // run the interface
    SP<IT> interface = new IT(arg);

  // set problem variables for this cycle
    delta_t = interface->get_delta_t();
  // dump_f  = interface->get_print_f();
    rnd_con = new Rnd_Control(interface->get_seed());
    Particle_Buffer<PT>::set_buffer_size(interface->get_buffer());

  // initialize the mesh builder and build mesh
    {
	BT mt_builder(interface);
	mesh = mt_builder.build_Mesh();
    }

  // initialize the opacity builder and build the state
    {
	Opacity_Builder<MT> opacity_builder(interface);

      // build the mat_state and opacity
	mat_state = opacity_builder.build_Mat(mesh);
	opacity   = opacity_builder.build_Opacity(mesh, mat_state);
    }

  // do the source initialization
    source_init = new Parallel_Source_Init<MT>(interface, mesh);
    if (cycle == 0)
    {
      // create the initial census
	IT::set_census(source_init->calc_initial_census
		       (mesh, opacity, mat_state, rnd_con)); 
    }
    else
    {
      // make a communications buffer
	buffer = new Particle_Buffer<PT>(*mesh, *rnd_con);

      // initialize the source
	source = source_init->initialize(mesh, opacity, mat_state, 
					 rnd_con, *buffer);
    
      // make a tally
	tally = new Tally<MT>(mesh, source_init->get_evol_net());
	
      // make sure tally, source, and buffer are made
	Ensure (tally);
	Ensure (buffer);
	Ensure (source);
	Ensure (IT::get_census()->size() == 0);
    }

  // reclaim what memory we can
    kill(source_init);

  // check to make sure we have all the necessary objects
    Ensure (mesh);
    Ensure (opacity);
    Ensure (mat_state);
    Ensure (rnd_con);
    Ensure (IT::get_census());
    Ensure (!source_init);    
}

//---------------------------------------------------------------------------//
// IMC TRANSPORT
//---------------------------------------------------------------------------//
// do a timestep on a FULLY REPLICATED MESH

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::step_IMC()
{
  // objects required for IMC transport
    Require (cycle > 0);
    Require (mesh);
    Require (opacity);
    Require (mat_state);
    Require (buffer);
    Require (source);
    Require (tally);

    cerr << ">> Doing IMC transport for cycle " << cycle
	 << " on proc " << node() << " using full replication." << endl;

  // make a new census Comm_Buffer on this node
    new_census_bank = new Particle_Buffer<PT>::Census();

  // diagnostic
    SP<typename PT::Diagnostic> check = new PT::Diagnostic(cout, true);

  // get source particles and run them to completion
    int counter = 0;
    while (*source)
    {
      // get a particle from the source
	SP<PT> particle = source->get_Source_Particle(delta_t);
	Check (particle->status());

      // transport the particle
	particle->transport(*mesh, *opacity, *tally);
	counter++;

      // after the particle is no longer active take appropriate action
	Check (!particle->status());
	
      // if census write to file
	if (particle->desc() == "census")
	    new_census_bank->push(particle);

      // message particle counter
	if (!(counter % dump_f)) 
	    cerr << setw(10) << counter << " particles run on proc " 
		 << node()   << endl;
    }

  // finished with this timestep
    cerr << ">> Finished particle transport for cycle " << cycle
	 << " on proc " << node() << endl;

  // now remove objects we no longer need
    kill (source);
    kill (mat_state);
    kill (opacity);
    kill (mesh);
    kill (rnd_con);
    kill (buffer);

  // object inventory
    Ensure (tally);	
    Ensure (new_census_bank);
    Ensure (!source);
    Ensure (!mat_state);
    Ensure (!opacity);
    Ensure (!source_init);
    Ensure (!mesh);
    Ensure (!rnd_con);
    Ensure (!buffer);
}

//---------------------------------------------------------------------------//
// ACCUMULATE TOTALS
//---------------------------------------------------------------------------//

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::regroup(typename IT::Arguments &arg)
{
    Require (tally);
    Require (new_census_bank);    

  // set new global census
    IT::set_census(new_census_bank);
    kill (new_census_bank);  

   // write tallies to return data
    Check (tally->num_cells() == arg.num_cells);
    for (int i = 1; i <= tally->num_cells(); i++)
    {
      	arg.e_dep[i-1]   = tally->get_energy_dep(i) - tally->get_evol_net(i); 
	arg.rad_den[i-1] = tally->get_accum_ewpl(i) / (Global::c *
						       tally->volume(i) *
						       delta_t);
    }

  // we are done
    Ensure (IT::get_census());
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Host_Manager.cc
//---------------------------------------------------------------------------//
