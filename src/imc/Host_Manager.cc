//----------------------------------*-C++-*----------------------------------//
// Host_Manager.cc
// Thomas M. Evans
// Sun Aug  2 12:43:05 1998
//---------------------------------------------------------------------------//
// @> Host_Manager class implementation file
//---------------------------------------------------------------------------//

#include "imc/Host_Manager.hh"
#include "imc/Parallel_Builder.hh"
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

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor

template<class MT, class BT, class IT, class PT>
Host_Manager<MT,BT,IT,PT>::Host_Manager()
    : delta_t(0), cycle(0), dump_f(Global::huge)
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
    Check (!Global::host_census);
}

//---------------------------------------------------------------------------//
// execute IMC in the host
//---------------------------------------------------------------------------//
// run an IMC cycle

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::execute_IMC(const typename IT::Arguments
					    &arg)
{
  // run through an IMC timestep
  
  // initialize IMC
    initialize(arg);

  // do a time-step
    step_IMC();

  // cleanup before shipping off
    regroup();
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
    dump_f  = interface->get_printf();
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
    source_init = new Source_Init<MT>(interface, mesh);
    source_init->set_census(Global::host_census);
    source_init->host_init(mesh, opacity, mat_state, rnd_con, cycle);
    if (!Global::host_census)
	Global::host_census = source_init->get_census();

  // make a communications buffer
    buffer = new Particle_Buffer<PT>(*mesh, *rnd_con);
    
  // make a parallel_builder (TEMPORARY)
    SP<Parallel_Builder<MT> > parallel_builder = new Parallel_Builder<MT>();
    source = parallel_builder->send_Source(mesh, mat_state, rnd_con,
					   *source_init, *buffer);
    Check (Global::host_census->size() == 0);

  // make a tally
    tally = new Tally<MT>(mesh);

  // reclaim what memory we can
    kill(source_init);
    
  // check to make sure we have all the necessary objects
    Ensure (mesh);
    Ensure (opacity);
    Ensure (mat_state);
    Ensure (rnd_con);
    Ensure (tally);
    Ensure (source);
    Ensure (buffer);
    Ensure (Global::host_census);
    Ensure (!source_init);
}

//---------------------------------------------------------------------------//
// do a timestep

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::step_IMC()
{}

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::regroup()
{}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Host_Manager.cc
//---------------------------------------------------------------------------//
