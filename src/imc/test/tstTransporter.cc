//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstTransporter.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:15:17 2001
 * \brief  Transporter classes test.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "tstTransporter.hh"
#include "../Transporter.hh"
#include "../Rep_Transporter.hh"
#include "../DD_Transporter.hh"
#include "../Gray_Particle.hh"
#include "../Multigroup_Particle.hh"
#include "../Frequency.hh"
#include "../Release.hh"
#include "../Mat_State.hh"
#include "../Opacity.hh"
#include "../Flat_Mat_State_Builder.hh"
#include "../Rep_Source_Builder.hh"
#include "../Source.hh"
#include "../Tally.hh"
#include "../Diffusion_Opacity.hh"
#include "../Random_Walk.hh"
#include "../Extrinsic_Surface_Tracker.hh"
#include "../Tally_Builder.hh"
#include "mc/Communicator.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Comm_Patterns.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <sstream>
#include <typeinfo>

using namespace std;

using rtt_imc_test::RW_Interface;
using rtt_imc_test::IMC_Flat_Interface;
using rtt_imc_test::Parser;
using rtt_imc::Transporter;
using rtt_imc::Rep_Transporter;
using rtt_imc::DD_Transporter;
using rtt_imc::Flat_Mat_State_Builder;
using rtt_imc::Source;
using rtt_imc::Opacity;
using rtt_imc::Mat_State;
using rtt_imc::Rep_Source_Builder;
using rtt_imc::Tally;
using rtt_imc::Tally_Builder;
using rtt_imc::Random_Walk;
using rtt_imc::Diffusion_Opacity;
using rtt_imc::Extrinsic_Surface_Tracker;
using rtt_mc::Communicator;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::Comm_Patterns;
using rtt_rng::Rnd_Control;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

// typedefs
typedef rtt_mc::OS_Mesh                       MT;
typedef rtt_imc::Gray_Frequency               Gray;
typedef rtt_imc::Multigroup_Frequency         MG;
typedef rtt_imc::Gray_Particle<OS_Mesh>       GPT;
typedef rtt_imc::Multigroup_Particle<OS_Mesh> MGPT; 

//---------------------------------------------------------------------------//
// SERVICES
//---------------------------------------------------------------------------//

template<class FT, class PT>
void check_tallies(SP<Tally<MT> > tally, int num_run)
{
    FAILMSG("We should not get this check.");
}

//---------------------------------------------------------------------------//

template<>
void check_tallies<Gray,GPT>(SP<Tally<MT> > tally, int num_run)
{
    rtt_c4::global_barrier();

    // do some integral checks on tally output
    double ew_escaped  = tally->get_ew_escaped();
    double erg_dep_tot = tally->get_energy_dep_tot();
    double ecen_tot    = tally->get_new_ecen_tot();
    int    ncen_tot    = tally->get_new_ncen_tot();
    int    neff_scat   = tally->get_accum_n_effscat();
    int    nbnd_cross  = tally->get_accum_n_bndcross();
    int    nescaped    = tally->get_accum_n_escaped();

    // do sums
    rtt_c4::global_sum(ew_escaped);
    rtt_c4::global_sum(erg_dep_tot);
    rtt_c4::global_sum(ecen_tot);
    rtt_c4::global_sum(ncen_tot);
    rtt_c4::global_sum(neff_scat);
    rtt_c4::global_sum(nbnd_cross);
    rtt_c4::global_sum(nescaped);

    // check sums
    if (!soft_equiv(erg_dep_tot, 2.293336, 1.e-6))  ITFAILS;
    if (!soft_equiv(ew_escaped, 107.165317, 1.e-6)) ITFAILS;
    if (!soft_equiv(ecen_tot, 1881.311183, 1.e-6))  ITFAILS;
    if (ncen_tot != 941)                            ITFAILS;
    if (neff_scat != 12)                            ITFAILS;
    if (nbnd_cross != 185)                          ITFAILS;
    if (nescaped != 54)                             ITFAILS;

    // check num_run
    if (num_run != 995) ITFAILS;
    
    if (rtt_imc_test::passed)
	PASSMSG("Integral transport checks ok for gray problem.")
}

//---------------------------------------------------------------------------//

template<>
void check_tallies<MG,MGPT>(SP<Tally<MT> > tally, int num_run)
{
    rtt_c4::global_barrier();

    // do some integral checks on tally output
    double ew_escaped  = tally->get_ew_escaped();
    double erg_dep_tot = tally->get_energy_dep_tot();
    double ecen_tot    = tally->get_new_ecen_tot();
    int    ncen_tot    = tally->get_new_ncen_tot();
    int    neff_scat   = tally->get_accum_n_effscat();
    int    nbnd_cross  = tally->get_accum_n_bndcross();
    int    nescaped    = tally->get_accum_n_escaped();

    // do sums
    rtt_c4::global_sum(ew_escaped);
    rtt_c4::global_sum(erg_dep_tot);
    rtt_c4::global_sum(ecen_tot);
    rtt_c4::global_sum(ncen_tot);
    rtt_c4::global_sum(neff_scat);
    rtt_c4::global_sum(nbnd_cross);
    rtt_c4::global_sum(nescaped);

    // check sums
    if (!soft_equiv(erg_dep_tot, 2.159577, 1.e-6))  ITFAILS;
    if (!soft_equiv(ew_escaped, 93.0673595, 1.e-6)) ITFAILS;
    if (!soft_equiv(ecen_tot, 1797.756481, 1.e-6))  ITFAILS;
    if (ncen_tot != 948)                            ITFAILS;
    if (neff_scat != 153)                           ITFAILS;
    if (nbnd_cross != 190)                          ITFAILS;
    if (nescaped != 49)                             ITFAILS;

    // check num_run
    if (num_run != 997) ITFAILS;
    
    if (rtt_imc_test::passed)
	PASSMSG("Integral transport checks ok for multigroup problem.")
}

//---------------------------------------------------------------------------//

template<class FT, class PT>
void check_rw_tallies(SP<Tally<MT> > tally, int num_run)
{
    FAILMSG("We should not get this check.");
}

//---------------------------------------------------------------------------//

template<>
void check_rw_tallies<Gray,GPT>(SP<Tally<MT> > tally, int num_run)
{
    rtt_c4::global_barrier();

    // do some integral checks on tally output
    double ew_escaped  = tally->get_ew_escaped();
    double erg_dep_tot = tally->get_energy_dep_tot();
    double ecen_tot    = tally->get_new_ecen_tot();
    int    ncen_tot    = tally->get_new_ncen_tot();
    int    neff_scat   = tally->get_accum_n_effscat();
    int    nbnd_cross  = tally->get_accum_n_bndcross();
    int    nescaped    = tally->get_accum_n_escaped();
    int    nrws        = tally->get_RW_Sub_Tally()->get_accum_n_random_walks();

    // do sums
    rtt_c4::global_sum(ew_escaped);
    rtt_c4::global_sum(erg_dep_tot);
    rtt_c4::global_sum(ecen_tot);
    rtt_c4::global_sum(ncen_tot);
    rtt_c4::global_sum(neff_scat);
    rtt_c4::global_sum(nbnd_cross);
    rtt_c4::global_sum(nescaped);
    rtt_c4::global_sum(nrws);

    // check sums
    if (!soft_equiv(erg_dep_tot, 0.91662919, 1.e-6)) ITFAILS;
    if (!soft_equiv(ew_escaped, 0.04447913, 1.e-6))  ITFAILS;
    if (!soft_equiv(ecen_tot, 1.7061015068, 1.e-6))  ITFAILS;
    if (ncen_tot != 980)                             ITFAILS;
    if (neff_scat != 1012)                           ITFAILS;
    if (nbnd_cross != 179)                           ITFAILS;
    if (nescaped != 20)                              ITFAILS;
    if (nrws   != 1301)                              ITFAILS;

    // check num_run
    if (num_run != 1000) ITFAILS;
    
    if (rtt_imc_test::passed)
	PASSMSG("Integral transport checks ok for gray random walk problem.")
}

//---------------------------------------------------------------------------//

template<>
void check_rw_tallies<MG,MGPT>(SP<Tally<MT> > tally, int num_run)
{
    rtt_c4::global_barrier();

    // do some integral checks on tally output
    double ew_escaped  = tally->get_ew_escaped();
    double erg_dep_tot = tally->get_energy_dep_tot();
    double ecen_tot    = tally->get_new_ecen_tot();
    int    ncen_tot    = tally->get_new_ncen_tot();
    int    neff_scat   = tally->get_accum_n_effscat();
    int    nbnd_cross  = tally->get_accum_n_bndcross();
    int    nescaped    = tally->get_accum_n_escaped();
    int    nrws        = tally->get_RW_Sub_Tally()->get_accum_n_random_walks();

    // do sums
    rtt_c4::global_sum(ew_escaped);
    rtt_c4::global_sum(erg_dep_tot);
    rtt_c4::global_sum(ecen_tot);
    rtt_c4::global_sum(ncen_tot);
    rtt_c4::global_sum(neff_scat);
    rtt_c4::global_sum(nbnd_cross);
    rtt_c4::global_sum(nescaped);
    rtt_c4::global_sum(nrws);

    // check sums
    if (!soft_equiv(erg_dep_tot, 0.081283934, 1.e-6)) ITFAILS;
    if (!soft_equiv(ew_escaped, 0.0080145985, 1.e-6)) ITFAILS;
    if (!soft_equiv(ecen_tot, 0.279047012801, 1.e-6)) ITFAILS;
    if (ncen_tot != 970)                              ITFAILS;
    if (neff_scat != 258)                             ITFAILS;
    if (nbnd_cross != 193)                            ITFAILS;
    if (nescaped != 26)                               ITFAILS;
    if (nrws   != 195)                                ITFAILS;

    // check num_run
    if (num_run != 996) ITFAILS;
    
    if (rtt_imc_test::passed)
	PASSMSG("Integral transport checks ok for mg random walk problem.")
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

template<class MT, class FT, class PT>
void rep_transporter_test()
{
    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<MT> mesh = mb->build_Mesh();

    // build a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a replication transporter
    SP<Transporter<MT,FT,PT> > transporter;
    transporter = new Rep_Transporter<MT,FT,PT>(topology);

    // should be full rep transporter
    if (transporter->type() != "replication") ITFAILS;

    // transporter should not be ready
    if (transporter->ready()) ITFAILS;
}

//---------------------------------------------------------------------------//

template<class MT, class FT, class PT>
void rep_transporter_run_test()
{
    // set the rn_stream to 0
    rtt_rng::rn_stream = 0;

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<MT> mesh = mb->build_Mesh();

    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));

    // build a topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // get an interface (dummy)
    SP<IMC_Flat_Interface<PT> > interface(new IMC_Flat_Interface<PT>(mb));

    // build the Mat_State
    Flat_Mat_State_Builder<MT,FT> ob(interface);
    ob.build_mat_classes(mesh);

    SP<FT>              frequency = ob.get_Frequency();
    SP<Mat_State<MT> >  mat       = ob.get_Mat_State();
    SP<Opacity<MT,FT> > opacity   = ob.get_Opacity();

    // build a Rep_Source Builder
    Rep_Source_Builder<MT,FT,PT> source_builder(interface, mesh, topology);

    // build the source
    SP<Source<MT,FT,PT> > source = 
	source_builder.build_Source(mesh, mat, opacity, rcon, patterns);

    // build a replication transporter
    SP<Transporter<MT,FT,PT> > transporter;
    transporter = new Rep_Transporter<MT,FT,PT>(topology);

    // build a tally
    SP<Tally<MT> > tally(new Tally<MT>(mesh));

    // build a "NULL" communicator
    SP<Communicator<PT> > comm;

    // build a random walk object
    SP<Random_Walk<MT> > rwalk;

    // build an extrinsic surface tracker
    SP<Extrinsic_Surface_Tracker> tracker;
    
    // set the transporter
    transporter->set(mesh, mat, opacity, source, tally, rwalk, tracker, comm);

    // transport
    double dt = interface->get_delta_t();
    int cycle = interface->get_cycle();
    transporter->transport(dt, cycle, 100000, 0, 0);

    // check tallies
    check_tallies<FT,PT>(tally, transporter->get_num_run());

    if (rtt_imc_test::passed)
    {
	ostringstream message;
	message << "Finished a replication transport cycle for " 
		<< typeid(FT).name() << " and " << typeid(PT).name()
		<< endl;
	PASSMSG(message.str().c_str());
    }
    else
    {
	ostringstream message;
	message << "Unable to complete a replication transport cycle for " 
		<< typeid(FT).name() << " and " << typeid(PT).name()
		<< endl;
	FAILMSG(message.str().c_str());
    }
}

//---------------------------------------------------------------------------//

template<class MT, class FT, class PT>
void rep_transporter_random_walk_run_test()
{
    // set the rn_stream to 0
    rtt_rng::rn_stream = 0;

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input_3D"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<MT> mesh = mb->build_Mesh();

    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));

    // build a topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // get an interface (dummy)
    SP<RW_Interface<PT> > interface(new RW_Interface<PT>(mb));

    // build the Mat_State
    Flat_Mat_State_Builder<MT,FT> ob(interface);
    ob.build_mat_classes(mesh);

    SP<FT>              frequency           = ob.get_Frequency();
    SP<Mat_State<MT> >  mat                 = ob.get_Mat_State();
    SP<Opacity<MT,FT> > opacity             = ob.get_Opacity();
    SP<Diffusion_Opacity<MT> > diff_opacity = ob.get_Diffusion_Opacity();

    // build a Rep_Source Builder
    Rep_Source_Builder<MT,FT,PT> source_builder(interface, mesh, topology);

    // build the source
    SP<Source<MT,FT,PT> > source = 
	source_builder.build_Source(mesh, mat, opacity, rcon, patterns);

    // build a replication transporter
    SP<Transporter<MT,FT,PT> > transporter;
    transporter = new Rep_Transporter<MT,FT,PT>(topology);

    // build a tally
    Tally_Builder<MT> tally_builder(interface);
    SP<Tally<MT> > tally = tally_builder.build_Tally(mesh);

    // build a "NULL" communicator
    SP<Communicator<PT> > comm;

    // build a random walk object
    SP<Random_Walk<MT> > rwalk(new Random_Walk<MT>(mesh, diff_opacity));
    if (!rwalk) ITFAILS;

    // build an extrinsic surface tracker
    SP<Extrinsic_Surface_Tracker> tracker;
    
    // set the transporter
    transporter->set(mesh, mat, opacity, source, tally, rwalk, tracker, comm);

    // transport
    double dt = interface->get_delta_t();
    int cycle = interface->get_cycle();
    transporter->transport(dt, cycle, 100000, 0, 0);

    // check tallies
    check_rw_tallies<FT,PT>(tally, transporter->get_num_run());

    if (rtt_imc_test::passed)
    {
	ostringstream message;
	message << "Finished replication transport cycle for " 
		<< typeid(FT).name() << " and " << typeid(PT).name()
		<< " that did " 
		<< tally->get_RW_Sub_Tally()->get_accum_n_random_walks()
		<< " random walks" << endl;
	PASSMSG(message.str().c_str());
    }
    else
    {
	ostringstream message;
	message << "Unable to complete a replication transport cycle for " 
		<< typeid(FT).name() << " and " << typeid(PT).name()
		<< endl;
	FAILMSG(message.str().c_str());
    }
}

//---------------------------------------------------------------------------//

template<class MT, class FT, class PT>
void DD_transporter_test()
{
    // only perform this test on two processors
    if (C4::nodes() != 2) 
	return;
    
    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<MT> mesh = mb->build_Mesh();

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

    // make a rnd control
    Rnd_Control rcon(347223);

    // build a DD transporter
    SP<Transporter<MT,FT,PT> > transporter;
    transporter = new DD_Transporter<MT,FT,PT>(topology);

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

	rep_transporter_test<MT,Gray,GPT>();
	DD_transporter_test<MT,Gray,GPT>();

	// run some particles
	rep_transporter_run_test<MT,Gray,GPT>();
	rep_transporter_run_test<MT,MG,MGPT>();

	// run some particles with random walk
	rep_transporter_random_walk_run_test<MT,Gray,GPT>();
	rep_transporter_random_walk_run_test<MT,MG,MGPT>();
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
