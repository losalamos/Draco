//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstTransporter.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:15:17 2001
 * \brief  Transporter classes test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Transporter.hh"
#include "../Rep_Transporter.hh"
#include "../DD_Transporter.hh"
#include "../Particle.hh"
#include "../Release.hh"
#include "../Mat_State.hh"
#include "../Opacity.hh"
#include "../Flat_Mat_State_Builder.hh"
#include "../Rep_Source_Builder.hh"
#include "../Source.hh"
#include "../Tally.hh"
#include "mc/Particle_Buffer.hh"
#include "mc/Communicator.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Parallel_Data_Operator.hh"
#include "mc/Comm_Patterns.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

using rtt_imc_test::IMC_Flat_Interface;
using rtt_imc_test::Parser;
using rtt_imc::Transporter;
using rtt_imc::Rep_Transporter;
using rtt_imc::DD_Transporter;
using rtt_imc::Particle;
using rtt_imc::Flat_Mat_State_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Opacity;
using rtt_imc::Source;
using rtt_imc::Rep_Source_Builder;
using rtt_imc::Tally;
using rtt_mc::Particle_Buffer;
using rtt_mc::Communicator;
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

//---------------------------------------------------------------------------//
// TESTS
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
    SP<IMC_Flat_Interface> interface(new IMC_Flat_Interface(mb));

    // build the Mat_State
    Flat_Mat_State_Builder<OS_Mesh> ob(interface);
    SP<OS_Mat> mat         = ob.build_Mat_State(mesh);
    SP<OS_Opacity> opacity = ob.build_Opacity(mesh, mat);

    // build a Rep_Source Builder
    Rep_Source_Builder<OS_Mesh> source_builder(interface, mesh, topology);

    // build the source
    SP<Source<OS_Mesh, Particle<OS_Mesh> > > source = 
	source_builder.build_Source(mesh, mat, opacity, rcon, patterns);

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
    SP<IMC_Flat_Interface> interface(new IMC_Flat_Interface(mb, 3));

    // build a DD topology
    SP<Topology> topology;
    if (C4::node() == 0)
    {
	if (mesh->num_cells() != 6) ITFAILS;
	
	// build a DD topology with cells 1,2,3 on proc 0 and cells 4,5,6 on
	// proc 1; we will expand this when we expand the test in general
	
	Topology::vf_int cpp(2, vector<int>(3));
	Topology::vf_int ppc(6, vector<int>(1));
	Topology::vf_int bc(2, vector<int>(3));

	// procs per cell
	ppc[0][0] = 0;
	ppc[1][0] = 0;
	ppc[2][0] = 0;
	ppc[3][0] = 1;
	ppc[4][0] = 1;
	ppc[5][0] = 1;

	// cells per proc
	cpp[0][0] = 1;
	cpp[0][1] = 2;
	cpp[0][2] = 3;
	cpp[1][0] = 4;
	cpp[1][1] = 5;
	cpp[1][2] = 6;

	// boundary cell data
	bc[0][0] = 4;
	bc[0][1] = 5;
	bc[0][2] = 6;
	bc[1][0] = 1;
	bc[1][1] = 2;
	bc[1][2] = 3;

	topology = new General_Topology(cpp, ppc, bc, "DD");

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
    Rnd_Control rcon(347223);

    // build a DD transporter
    SP<Transporter<MT,PT> > transporter;
    transporter = new DD_Transporter<MT,PT>(topology);

    // should be full DD transporter
    if (transporter->type() != "DD") ITFAILS;

    // transporter should not be ready
    if (transporter->ready()) ITFAILS;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_imc::release() 
		 << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	rep_transporter_test<OS_Mesh,Particle<OS_Mesh> >();
	DD_transporter_test<OS_Mesh, Particle<OS_Mesh> >();

	// run some particles
	rep_transporter_run_test<OS_Mesh,Particle<OS_Mesh> >();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstTransporter, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstTransporter Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstTransporter on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstTransporter.cc
//---------------------------------------------------------------------------//
