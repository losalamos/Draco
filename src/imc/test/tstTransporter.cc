//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstTransporter.cc
 * \author Thomas M. Evans
 * \date   Mon Apr 17 11:29:58 2000
 * \brief  Transporter class testing file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Transporter.hh"
#include "../Rep_Transporter.hh"
#include "../DD_Transporter.hh"
#include "../Particle_Buffer.hh"
#include "../Particle.hh"
#include "../Topology_Builder.hh"
#include "../Release.hh"
#include "../Mat_State.hh"
#include "../Opacity.hh"
#include "../Opacity_Builder.hh"
#include "../Rep_Source_Builder.hh"
#include "../Source.hh"
#include "../Tally.hh"
#include "../Communicator.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Parallel_Data_Operator.hh"
#include "mc/Comm_Patterns.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

using rtt_imc_test::IMC_Interface;
using rtt_imc_test::Parser;
using rtt_imc::Transporter;
using rtt_imc::Rep_Transporter;
using rtt_imc::DD_Transporter;
using rtt_imc::Particle_Buffer;
using rtt_imc::Particle;
using rtt_imc::Topology_Builder;
using rtt_imc::Opacity_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Opacity;
using rtt_imc::Source;
using rtt_imc::Rep_Source_Builder;
using rtt_imc::Communicator;
using rtt_imc::Tally;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::Parallel_Data_Operator;
using rtt_mc::Comm_Patterns;
using rtt_rng::Rnd_Control;
using rtt_dsxx::SP;

// typedefs
typedef Mat_State<OS_Mesh> OS_Mat;
typedef Opacity<OS_Mesh>   OS_Opacity;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

//---------------------------------------------------------------------------//

template<class MT, class PT>
void rep_transporter_test()
{
    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a replication transporter
    SP<Transporter<MT,PT> > transporter;
    transporter = new Rep_Transporter<MT,PT>(topology);

    // should be full rep transporter
    if (transporter->type() != "replication") ITFAILS;

    // transporter should not be ready
    if (transporter->ready()) ITFAILS;
}

//---------------------------------------------------------------------------//

template<class MT, class PT>
void rep_transporter_run_test()
{
    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));

    // build a topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // get an interface (dummy)
    SP<IMC_Interface> interface(new IMC_Interface(mb));

    // build the Mat_State
    Opacity_Builder<OS_Mesh> ob(interface);
    SP<OS_Mat> mat         = ob.build_Mat(mesh);
    SP<OS_Opacity> opacity = ob.build_Opacity(mesh, mat);

    // build a Rep_Source Builder
    Rep_Source_Builder<OS_Mesh> source_builder(interface, mesh, topology);

    // build the source
    SP<Source<OS_Mesh> > source = source_builder.build_Source(mesh, mat,
							      opacity, rcon,
							      patterns);

    cout << source->get_num_source_particles() << endl;

    // build a replication transporter
    SP<Transporter<MT,PT> > transporter;
    transporter = new Rep_Transporter<MT,PT>(topology);

    // build a tally
    SP<Tally<OS_Mesh> > tally(new Tally<OS_Mesh>(mesh));

    // build a "NULL" communicator
    SP<Communicator<PT> > comm;
    
    // set the transporter
    transporter->set(mesh, mat, opacity, source, tally, comm);

    // transport
    double dt = interface->get_delta_t();
    int cycle = interface->get_cycle();
    transporter->transport(dt, cycle, 100000, 0, 0);
}

//---------------------------------------------------------------------------//

template<class MT, class PT>
void DD_transporter_test()
{
    // only perform this test on two processors
    if (C4::nodes() != 2) 
	return;
    
    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // get the dummy interface with a capacity of 3 cells (2 processor)
    SP<IMC_Interface> interface(new IMC_Interface(mb, 3));

    // build the Topology builder and full replication topology
    SP<Topology> topology;
    if (C4::node() == 0)
    {
	Topology_Builder<OS_Mesh> tb(interface);
	topology = tb.build_Topology(mesh);

	// cast to general topology for sending
	const General_Topology *gt = dynamic_cast
	    <General_Topology *>(topology.bp());
	if (!gt) ITFAILS;

	rtt_imc_test::send_TOP(*gt);
    }
    if (C4::node())
    {
	// receive the topology
	topology = rtt_imc_test::recv_TOP();
    }

    // check to make sure Topology is full DD
    if (topology->get_parallel_scheme() != "DD") ITFAILS;

    // make a Particle_Buffer
    Rnd_Control rcon(324235);
    SP<Particle_Buffer<PT> > buffer(new Particle_Buffer<PT>(*mesh, rcon));

    // build a DD transporter
    SP<Transporter<MT,PT> > transporter;
    transporter = new DD_Transporter<MT,PT>(topology, buffer);

    // should be full DD transporter
    if (transporter->type() != "DD") ITFAILS;

    // transporter should not be ready
    if (transporter->ready()) ITFAILS;
}

//---------------------------------------------------------------------------//
// main

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl; 
	    C4::Finalize();
	    return 0;
	}

    // tests
    try
    {
	rep_transporter_test<OS_Mesh,Particle<OS_Mesh> >();
	DD_transporter_test<OS_Mesh, Particle<OS_Mesh> >();

	// run some particles
	rep_transporter_run_test<OS_Mesh,Particle<OS_Mesh> >();
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	cout << "Test: assertion failure at line " 
	     << ass.what() << endl;
	C4::Finalize();
	return 1;
    }
    catch(...)
    {
	cout << "HELP ME" << endl;
	C4::Finalize();
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "***************************************" << endl; 
    if (passed) 
    {
        cout << "**** Transporter Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "***************************************" << endl;
    cout << endl;

    cout << "Done testing Transporter on node: " << C4::node() << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstTransporter.cc
//---------------------------------------------------------------------------//
