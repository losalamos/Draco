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

IMCSPACE

// draco components
using C4::node;
using C4::nodes;
using C4::HTSyncSpinLock;

// stl components
using std::cout;
using std::cerr;
using std::endl;

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
	cout << "** Read input file " << argv << " on node " << node()
	     << endl;

  // make the rnd_control object and set the buffer sizes
    rnd_con = new Rnd_Control(9836592);
    Particle_Buffer<PT>::set_buffer_size(100);

  // initialize the mesh builder and build mesh
    {
	BT mt_builder(interface);
	mesh = mt_builder.build_Mesh();
	if (verbose)
	    cout << "** Built mesh on node " << node() << endl;
    }

  // initialize the opacity builder and build the state
    {
	Opacity_Builder<MT>  opacity_builder(interface, mesh);
	mat_state    = opacity_builder.build_Mat();
	global_state = new Mat_State<MT>::Shell(*mat_state);
	opacity      = opacity_builder.build_Opacity();
	if (verbose)
	    cout << "** Built Mat_State and Opacity on node " << node() 
		 << endl; 
    }

  // do the source initialization
    source_init = new Source_Init<MT>(interface, mesh);
    source_init->initialize(mesh, opacity, mat_state, rnd_con, cycle);
    if (verbose)
	cout << "** Did the source initialization on node " << node()
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
      // make the Parallel Builder and Particle Buffer
	parallel_builder = new Parallel_Builder<MT>(*mesh, *source_init);
	buffer           = new Particle_Buffer<PT>(*mesh, *rnd_con);

      // send the buffer sizes and Rnd seed
	for (int i = 1; i < nodes(); i++)
	{
	    Send (rnd_con->get_seed(), i, 200);
	    Send (Particle_Buffer<PT>::get_buffer_s(), i, 201);
	    Send (buffer->get_dsize(), i, 202);	
	    Send (buffer->get_isize(), i, 203);
	    Send (buffer->get_csize(), i, 204);
	}

      // now send out the mesh
    }

  // lets receive stuff on the IMC nodes
    if (node())
    {
      // make the Rnd_Control, Parallel_Builder, and Particle_Buffer on this
      // node 
	{
	  // receive the data from the host
	    int seed;
	    int s, d, i, c;
	    Recv (seed, 0, 200);
	    Recv (s, 0, 201);
	    Recv (d, 0, 202);
	    Recv (i, 0, 203);
	    Recv (c, 0, 204);

	  // now make the objects
	    rnd_con = new Rnd_Control(seed);
	    Particle_Buffer<PT>::set_buffer_size(s);
	    parallel_builder = new Parallel_Builder<MT>();
	    buffer           = new Particle_Buffer<PT>(d, i, c);
	}
    }

    HTSyncSpinLock h;
    {
	cout << "Node " << node() << endl;
	cout << "Rnd_Seed " << rnd_con->get_seed() << endl;
	cout << "Rnd_Number " << rnd_con->get_number() << endl;
	cout << "Rnd_Stream " << rnd_con->get_num() << endl;
	cout << "Buffer size " << Particle_Buffer<PT>::get_buffer_s() <<
	    endl;
	cout << "Double size " << Particle_Buffer<PT>::get_buffer_d() <<
	    endl;
	cout << "Int size " << Particle_Buffer<PT>::get_buffer_i() <<
	    endl;
    	cout << "Char size " << Particle_Buffer<PT>::get_buffer_c() <<
	    endl;
	RNG::rn_stream = node();
	cout << "RN_Stream " << RNG::rn_stream << endl << endl;
    }
}


CSPACE

//---------------------------------------------------------------------------//
//                              end of IMC_Man.cc
//---------------------------------------------------------------------------//
