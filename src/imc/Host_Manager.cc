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
void Host_Manager<MT,BT,IT,PT>::execute_IMC(const typename IT::Arguments
					    &arg)
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

    ofstream output("mesh.out");
  // initialize the mesh builder and build mesh
    {
	BT mt_builder(interface);
	mesh = mt_builder.build_Mesh();
	output << *mesh << endl;
    }

  // initialize the opacity builder and build the state
    {
	Opacity_Builder<MT> opacity_builder(interface);

      // build the mat_state and opacity
	mat_state = opacity_builder.build_Mat(mesh);
	output << *mat_state << endl;
	opacity   = opacity_builder.build_Opacity(mesh, mat_state);
	output << *opacity << endl;
    }

  // do the source initialization
    source_init = new Source_Init<MT>(interface, mesh);
    if (IT::get_census())
	source_init->set_census(IT::get_census());
    source_init->host_init(mesh, opacity, mat_state, rnd_con, cycle);
    if (!IT::get_census())
	IT::set_census(source_init->get_census());

  // NEED TO FIX CENSUS ENERGY UPDATES, READ THROUGH CENSUS EACH CYCLE AND
  // FIX THIS UP

    cout << "The number of particles is " <<
	source_init->get_ncentot()+source_init->get_nvoltot()+
	source_init->get_nsstot() << endl;
    cout << "Census " << source_init->get_ncentot() << "\t" <<
	source_init->get_ecentot() << endl;
    cout << "Volume " << source_init->get_nvoltot() << "\t" << 
	source_init->get_evoltot() << endl;
    cout << "SS     " << source_init->get_nsstot() << "\t" <<
	source_init->get_esstot() << endl;
    cout << "The random stream is at: " << RNG::rn_stream << endl;

  // return if this is an initialization cycle
    if (cycle == 0)
	return;

  // make a communications buffer
    buffer = new Particle_Buffer<PT>(*mesh, *rnd_con);
    
  // make a parallel_builder (TEMPORARY)
    SP<Parallel_Builder<MT> > parallel_builder =
	new Parallel_Builder<MT>(*mesh, *source_init);
    Check (parallel_builder->get_parallel_scheme() == "replication");
    source = parallel_builder->send_Source(mesh, mat_state, rnd_con,
					   *source_init, *buffer);
    Check (IT::get_census()->size() == 0);

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
    Ensure (IT::get_census());
    Ensure (!source_init);
}

//---------------------------------------------------------------------------//
// do a timestep

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::step_IMC()
{
    cout << "Looks like we made it!" << endl;
    cout << "The number of particles is " <<
	source->get_ncentot()+source->get_nvoltot()+source->get_nsstot() <<
	endl;
    cout << "Census " << source->get_ncentot() << endl;
    cout << "Volume " << source->get_nvoltot() << endl;
    cout << "SS     " << source->get_nsstot() << endl;
    cout << "The random stream is at: " << RNG::rn_stream << endl;
}

template<class MT, class BT, class IT, class PT>
void Host_Manager<MT,BT,IT,PT>::regroup()
{}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Host_Manager.cc
//---------------------------------------------------------------------------//
