//----------------------------------*-C++-*----------------------------------//
// IMC_Man.cc
// Thomas M. Evans
// Wed Jun  3 10:36:11 1998
//---------------------------------------------------------------------------//
// @> IMC_Man class implementation file.
//---------------------------------------------------------------------------//

#include "imc/IMC_Man.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

IMCSPACE

// draco components
using C4::node;
using C4::nodes;
using C4::HTSyncSpinLock;

// stl components
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor

template<class MT, class BT, class IT, class PT>
IMC_Man<MT,BT,IT,PT>::IMC_Man(bool verbose_)
    : delta_t(0), cycle(1), max_cycle(0), verbose(verbose_)
{
  // all the SPs should be null defined
    Check (!mesh);
    Check (!opacity);
    Check (!mat_state);
    Check (!rnd_con);
    Check (!parallel_builder);
    Check (!buffer);
    Check (!source_init);
    Check (!global_state);
}

//---------------------------------------------------------------------------//
// host processor initialization
//---------------------------------------------------------------------------//
// here we build the Mesh, Opacity, Mat_State, and Source_Init on the host
// processor 

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::host_init(char *argv)
{
  // run only on the host node
    if (node())
	return;

  // build the objects on the host processor
    if (verbose)
	cout << ">> Running through problem builders on proc " << node()
	     << endl;
    
  // run the interface (parser to input)
    string infile = argv;
    SP<IT> interface = new IT(infile);
    interface->parser();
    if (verbose)
	cout << " ** Read input file " << argv << " on node " << node()
	     << endl;

  // make the rnd_control object, set the buffer size, and assign delta_t on
  // the first cycle
    if (!rnd_con)
    {
	delta_t = interface->get_delta_t();
	rnd_con = new Rnd_Control(9836592);
	Particle_Buffer<PT>::set_buffer_size(100);
    }

  // initialize the mesh builder and build mesh
    {
	BT mt_builder(interface);
	mesh = mt_builder.build_Mesh();
	if (verbose)
	    cout << " ** Built mesh on node " << node() << endl;
    }

  // initialize the opacity builder and build the state
    {
	Opacity_Builder<MT>  opacity_builder(interface);
	mat_state    = opacity_builder.build_Mat(mesh);
	opacity      = opacity_builder.build_Opacity(mesh, mat_state);	
	if (verbose)
	    cout << " ** Built Mat_State and Opacity on node " << node() 
		 << endl; 

      // build the global state if we are on the first cycle
	if (!global_state)
	    global_state = new Mat_State<MT>::Shell(*mat_state);
    }

  // do the source initialization
    source_init = new Source_Init<MT>(interface, mesh);
    source_init->initialize(mesh, opacity, mat_state, rnd_con, cycle);
    if (verbose)
	cout << " ** Did the source initialization on node " << node()
	     << endl << endl;
}

//---------------------------------------------------------------------------//
// IMC processor initialization
//---------------------------------------------------------------------------//
// here we make communication objects and send stuff out to the IMC problem
// processors

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::IMC_init()
{
  // first lets do stuff on the host node
    if (!node())
    {
      // make sure we have performed the host_init
	Require (mesh);
	Require (opacity);
	Require (mat_state);
	Require (global_state);
	Require (source_init);

      // make the Parallel Builder and Particle Buffer on the first cycle
	if (!parallel_builder)
	{
	    parallel_builder = new Parallel_Builder<MT>(*mesh, *source_init); 
	    buffer           = new Particle_Buffer<PT>(*mesh, *rnd_con);

	  // send the buffer sizes, Rnd seed, and delta_t
	    for (int i = 1; i < nodes(); i++)
	    {
		Send (rnd_con->get_seed(), i, 200);
		Send (Particle_Buffer<PT>::get_buffer_s(), i, 201);
		Send (buffer->get_dsize(), i, 202);	
		Send (buffer->get_isize(), i, 203);
		Send (buffer->get_csize(), i, 204);
		Send (delta_t, i, 205);
	    }
	}

      // now send out the mesh, Opacity (full rep right now)
	parallel_builder->send_Mesh(*mesh);
	parallel_builder->send_Opacity(*opacity);

      // send out the Mat_State
	mat_state = parallel_builder->send_Mat(mesh, *mat_state);

      // send out the Source
	source = parallel_builder->send_Source(mesh, mat_state, rnd_con,
					       *source_init, *buffer);

      // kill objects we no longer need
	kill(source_init);
    }

  // lets receive stuff on the IMC nodes
    if (node())
    {
      // make the Rnd_Control, Parallel_Builder, and Particle_Buffer on this
      // processor if we are on the first cycle
	if (!parallel_builder)
	{
	  // receive the data from the host proc the first cycle
	    int seed;
	    int s, d, i, c;
	    Recv (seed, 0, 200);
	    Recv (s, 0, 201);
	    Recv (d, 0, 202);
	    Recv (i, 0, 203);
	    Recv (c, 0, 204);
	    Recv (delta_t, 0, 205);
	    
	  // now make the objects
	    rnd_con = new Rnd_Control(seed);
	    Particle_Buffer<PT>::set_buffer_size(s);
	    parallel_builder = new Parallel_Builder<MT>();
	    buffer           = new Particle_Buffer<PT>(d, i, c);
	}
	
      // receive the objects necessary for transport
	mesh      = parallel_builder->recv_Mesh();
	opacity   = parallel_builder->recv_Opacity(mesh);
	mat_state = parallel_builder->recv_Mat(mesh);
	source    = parallel_builder->recv_Source(mesh, mat_state, rnd_con, 
						  *buffer);
	if (verbose)
	    cout << " ** We received the mesh, opacity, mat_state," 
		 << " and source on node " << node() << endl;
    }
    
  // make sure each processor has the requisite objects to do transport
    Ensure (mesh);
    Ensure (opacity);
    Ensure (mat_state);
    Ensure (source);
    Ensure (buffer);
    Ensure (parallel_builder);
    Ensure (rnd_con);
    Ensure (!source_init);
    Ensure (delta_t > 0);

  // print a message indicating success
    if (verbose)
	cout << " ** Node " << node() << " has all objects necessary" 
	     << " for IMC transport." << endl;
}

//---------------------------------------------------------------------------//
// run particles on the problem geometry

template<class MT, class BT, class IT, class PT>
void IMC_Man<MT,BT,IT,PT>::step_IMC()
{
    if (verbose)
	cout << ">> Doing particle transport for cycle " << cycle
	     << " on proc " << node() << endl;

  // make a new census file on this node
    ostringstream cen_title;
    cen_title << "new_census." << node();
    string cen_file = cen_title.str();
    ofstream census(cen_file.c_str());

  // make a tally object and diagnostic
    Tally<MT> tally(mesh);
    SP<typename PT::Diagnostic> check = new PT::Diagnostic(cout, true);
    cout << *source << endl;

  // get source particles and run them to completion
    while (*source)
    {
      // get a particle from the source
	SP<PT> particle = source->get_Source_Particle(delta_t);
	Check (particle->status());

      // transport the particle
	particle->transport(*mesh, *opacity, tally, check);

      // after the particle is no longer active take appropriate action
	Check (!particle->status());
	
      // if census write to file
	if (particle->desc() == "census")
	    buffer->write_census(census, *particle);
      // no other choices right now
    }

    HTSyncSpinLock h;
    {
	cout << ">> Tally on " << node() << endl;
	cout << tally << endl;
    }

  // <<CONTINUE HERE>>
  // 1) in-cycle census
  // 2) concatening the tallies
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of IMC_Man.cc
//---------------------------------------------------------------------------//
