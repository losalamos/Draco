//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstSource_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Dec  8 16:39:33 1999
 * \brief  Source_Builder test file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Source_Builder.hh"
#include "../Rep_Source_Builder.hh"
#include "../Opacity_Builder.hh"
#include "../Opacity.hh"
#include "../Mat_State.hh"
#include "../Topology_Builder.hh"
#include "../Source.hh"
#include "../Release.hh"
#include "../Particle.hh"
#include "../Source_Init.hh"
#include "../Global.hh"
#include "mc/Rep_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Parallel_Data_Operator.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

using rtt_imc_test::IMC_Interface;
using rtt_imc_test::Parser;
using rtt_imc::Source_Builder;
using rtt_imc::Rep_Source_Builder;
using rtt_imc::Opacity;
using rtt_imc::Opacity_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Source;
using rtt_imc::Topology_Builder;
using rtt_imc::Particle;
using rtt_imc::Source_Init;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::Parallel_Data_Operator;
using rtt_rng::Rnd_Control;
using rtt_dsxx::SP;

typedef SP<Particle<OS_Mesh> > SP_Particle;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

//---------------------------------------------------------------------------//
// source init test --> satisfies all topologies because source init must be
// run with a full mesh

void source_init_test()
{
    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build an interface to a six cell fully replicated mesh
    SP<IMC_Interface> interface(new IMC_Interface(mb));  

    // build a Source_Init object->doesn't do anything yet
    Source_Init<IMC_Interface, OS_Mesh> source_init(interface);
}

//---------------------------------------------------------------------------//
// build source test for a full replication topology --> tests
// Rep_Source_Builder 

void source_replication_test()
{
    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build an interface to a six cell fully replicated mesh
    SP<IMC_Interface> interface(new IMC_Interface(mb));

    // build a Source_Init object; however, we haven't added this
    // functionality yet so we don't do anything else with it
    Source_Init<IMC_Interface, OS_Mesh> source_init(interface);

    // build a Topology: we do not use the Topology builder here because the
    // topology builder is designed to work on the host processor only -->
    // instead we will just build a Replication topology on each mesh
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));
    Parallel_Data_Operator pop(topology);

    // build a Mat_State and Opacity
    Opacity_Builder<OS_Mesh> ob(interface);
    SP<Mat_State<OS_Mesh> > mat    = ob.build_Mat(mesh);
    SP<Opacity<OS_Mesh> > opacity  = ob.build_Opacity(mesh, mat);

    // build a Rep_Source Builder
    Rep_Source_Builder<OS_Mesh> source_builder(interface, mesh, topology);

    // build the source
    SP<Source<OS_Mesh> > source = source_builder.build_Source(mesh, mat,
							      opacity, rcon);

    // get the global numbers for each species
    int global_nsstot  = source_builder.get_nsstot();
    int global_ncentot = source_builder.get_ncentot();
    int global_nvoltot = source_builder.get_nvoltot();

    // check to make sure globals and locals match
    {
	int local_nsstot  = source->get_nsstot();
	int local_ncentot = source->get_ncentot();
	int local_nvoltot = source->get_nvoltot();
	C4::gsum(local_nsstot);
	C4::gsum(local_ncentot);
	C4::gsum(local_nvoltot);

	if (local_nsstot != global_nsstot)   ITFAILS;
	if (local_ncentot != global_ncentot) ITFAILS;
	if (local_nvoltot != global_nvoltot) ITFAILS;
    }

    // calculate random number arrays on each processor
    vector<int> rn_vol(mesh->num_cells());
    vector<int> rn_cen(mesh->num_cells());
    vector<int> rn_ss(mesh->num_cells());

    // numbers of particles per cell per species; NOTE: the census numbers
    // are from the initial census calculation and are used to get the
    // initial random number stream ids, these are not the final census
    // numbers
    vector<int> global_nvol(mesh->num_cells(), 1);
    vector<int> global_ncen(mesh->num_cells(), 137);
    vector<int> global_nss(mesh->num_cells(), 0);
    {
	global_nvol[3] = 2;
	global_nvol[4] = 2;
	global_nvol[5] = 2;
	global_nss[0]  = 82;
	global_nss[1]  = 82;
    }

    // local source numbers per processor per species
    vector<int> local_nvol(mesh->num_cells(), 0);
    vector<int> local_ncen(mesh->num_cells(), 0);
    vector<int> local_nss(mesh->num_cells(), 0);

    // sum of ids for surface source and volume source per cell
    vector<int> id_sum_ss(mesh->num_cells(), 0);
    vector<int> id_sum_vol(mesh->num_cells(), 0);

    int running_rn  = 0;
    int rn_marker   = 0;

    // calculate census random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread = global_ncen[cell-1] / C4::nodes();
	int leftover    = global_ncen[cell-1] - even_spread * C4::nodes();

	rn_cen[cell-1] = running_rn + even_spread * C4::node() +
	    rtt_mc::global::min(C4::node(), leftover);

	running_rn += global_ncen[cell-1];

	// pre-combed, pre total source iteration census numbers
	local_ncen[cell-1] = even_spread;

	if (C4::node() < leftover)
	    local_ncen[cell-1]++;
    }

    // set random number id marker to first volume source id
    rn_marker = running_rn;

    // calculate volume random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread = global_nvol[cell-1] / C4::nodes();
	int leftover    = global_nvol[cell-1] - even_spread * C4::nodes();

	rn_vol[cell-1] = running_rn + even_spread * C4::node() +
	    rtt_mc::global::min(C4::node(), leftover);

	running_rn += global_nvol[cell-1];

	local_nvol[cell-1] = even_spread;

	if (C4::node() < leftover)
	    local_nvol[cell-1]++;
    }

    // add up random number ids in each cell for volume source
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	for (int np = 0; np < global_nvol[cell-1]; np++)
	    id_sum_vol[cell-1] += rn_marker + np;

	rn_marker += global_nvol[cell-1];
    }

    Check (rn_marker == running_rn);

    // calculate surface source random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread = global_nss[cell-1] / C4::nodes();
	int leftover    = global_nss[cell-1] - even_spread * C4::nodes();

	rn_ss[cell-1] = running_rn + even_spread * C4::node() +
	    rtt_mc::global::min(C4::node(), leftover);

	running_rn += global_nss[cell-1];

	local_nss[cell-1] = even_spread;

	if (C4::node() < leftover)
	    local_nss[cell-1]++;
    }

    // add up random number ids in each cell for surface source
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	for (int np = 0; np < global_nss[cell-1]; np++)
	    id_sum_ss[cell-1] += rn_marker + np;

	rn_marker += global_nss[cell-1];
    }

    // get particles out of the source and add up the energy weights
    double cen_ew = 0.0;
    double vol_ew = 0.0;
    double ss_ew  = 0.0;

    // calculated sums of surface source and volume ids
    vector<int> calc_vol_rn_sum(mesh->num_cells(), 0);
    vector<int> calc_ss_rn_sum(mesh->num_cells(), 0);
    
    // get surface sources
    for (int i = 0; i < source->get_nsstot(); i++)
    {
	SP_Particle particle = source->get_Source_Particle(.001);
	ss_ew               += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (rn_ss[cell-1] > streamid)                      ITFAILS;
	if (streamid >= rn_ss[cell-1] + local_nss[cell-1]) ITFAILS;

	// add up ids
	calc_ss_rn_sum[cell-1] += streamid;
    }

    // get volume sources
    for (int i = 0; i < source->get_nvoltot(); i++)
    {
	SP_Particle particle = source->get_Source_Particle(.001);
	vol_ew              += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (rn_vol[cell-1] > streamid)                       ITFAILS;
	if (streamid >= rn_vol[cell-1] + local_nvol[cell-1]) ITFAILS;

	// add up ids
	calc_vol_rn_sum[cell-1] += streamid;
    }

    // sum up ids from each processor
    pop.global_sum(calc_ss_rn_sum.begin(), calc_ss_rn_sum.end());
    pop.global_sum(calc_vol_rn_sum.begin(), calc_vol_rn_sum.end());

    if (calc_ss_rn_sum != id_sum_ss)   ITFAILS;
    if (calc_vol_rn_sum != id_sum_vol) ITFAILS;
   
    // get census sources
    for (int i = 0; i < source->get_ncentot(); i++)
    {
	SP_Particle particle = source->get_Source_Particle(.001);
	cen_ew              += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (rn_cen[cell-1] > streamid)                       ITFAILS;
	if (streamid >= rn_cen[cell-1] + local_ncen[cell-1]) ITFAILS;
    }

    // the source should be empty
    if (*source) ITFAILS;

    // sum up sampled energy weights across processors
    C4::gsum(ss_ew);
    C4::gsum(vol_ew);
    C4::gsum(cen_ew);

    // subtract energy losses
    ss_ew  -= source_builder.get_eloss_ss();
    vol_ew -= source_builder.get_eloss_vol();
    cen_ew -= source_builder.get_eloss_cen();

    // compare to totals in source_builder
    
    // get deterministically calculated global values of total energy from
    // the source builder; NOTE: this is the initial census energy
    double ref_vol = source_builder.get_evoltot();
    double ref_ss  = source_builder.get_esstot();
    double ref_cen = source_builder.get_initial_census_energy();

    if (fabs(vol_ew - ref_vol) > 1.e-8 * ref_vol) ITFAILS;
    if (fabs(ss_ew - ref_ss) > 1.e-8 * ref_ss)    ITFAILS;
    if (fabs(cen_ew - ref_cen) > 1.e-8 * ref_cen) ITFAILS;

    // check to hand calculations of same energies
    double hand_vol = 15.317634;
    double hand_ss  = 329.052202;
    double hand_cen = 1646.400000;

    if (fabs(vol_ew - hand_vol) > 1.e-4 * hand_vol) ITFAILS;
    if (fabs(ss_ew - hand_ss) > 1.e-4 * hand_ss)    ITFAILS;
    if (fabs(cen_ew - hand_cen) > 1.e-4 * hand_cen) ITFAILS;

    // check mat_vol_src
    double mvs;

    mvs = source_builder.get_mat_vol_src(1);
    if (fabs(mvs - 0.011460) > 1.e-4 * 0.011460) ITFAILS;
    mvs = source_builder.get_mat_vol_src(2);
    if (fabs(mvs - 0.011460) > 1.e-4 * 0.011460) ITFAILS;
    mvs = source_builder.get_mat_vol_src(3);
    if (fabs(mvs - 0.011460) > 1.e-4 * 0.011460) ITFAILS;

    mvs = source_builder.get_mat_vol_src(4);
    if (fabs(mvs - 0.026382) > 1.e-4 * 0.026382) ITFAILS;
    mvs = source_builder.get_mat_vol_src(5);
    if (fabs(mvs - 0.026382) > 1.e-4 * 0.026382) ITFAILS;
    mvs = source_builder.get_mat_vol_src(6);
    if (fabs(mvs - 0.026382) > 1.e-4 * 0.026382) ITFAILS;

    mvs = source_builder.get_mat_vol_srctot();
    if (fabs(mvs - 0.113524) > 1.e-4 * .113524) ITFAILS;

    // check evol_net
    double evn;

    evn = source_builder.get_evol_net(1);
    if (fabs(evn - 0.471351) > 1.e-4 * 0.471351) ITFAILS;
    evn = source_builder.get_evol_net(2);
    if (fabs(evn - 0.471351) > 1.e-4 * 0.471351) ITFAILS;
    evn = source_builder.get_evol_net(3);
    if (fabs(evn - 0.471351) > 1.e-4 * 0.471351) ITFAILS;

    evn = source_builder.get_evol_net(4);
    if (fabs(evn - 3.472368) > 1.e-4 * 3.472368) ITFAILS;
    evn = source_builder.get_evol_net(5);
    if (fabs(evn - 3.472368) > 1.e-4 * 3.472368) ITFAILS;
    evn = source_builder.get_evol_net(6);
    if (fabs(evn - 3.472368) > 1.e-4 * 3.472368) ITFAILS;
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

    // check source init -> independent of topology
    source_init_test();

    // full replication source test
    source_replication_test();

    // status of test
    cout << endl;
    cout <<     "******************************************" << endl; 
    if (passed) 
    {
        cout << "**** Source_Builder Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "******************************************" << endl;
    cout << endl;

    cout << "Done testing Source_Builder on node: " << C4::node() << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstSource_Builder.cc
//---------------------------------------------------------------------------//
