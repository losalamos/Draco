//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstSource_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:09:04 2001
 * \brief  Source_Builder and Source tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "IMC_DD_Test.hh"
#include "../Source_Builder.hh"
#include "../Rep_Source_Builder.hh"
#include "../DD_Source_Builder.hh"
#include "../Flat_Mat_State_Builder.hh"
#include "../Opacity.hh"
#include "../Mat_State.hh"
#include "../Source.hh"
#include "../Release.hh"
#include "../Frequency.hh"
#include "../Gray_Particle.hh"
#include "../Multigroup_Particle.hh"
#include "../Global.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Parallel_Data_Operator.hh"
#include "mc/Math.hh"
#include "mc/Comm_Patterns.hh"
#include "cdi/CDI.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/SP.hh"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

using rtt_imc_test::IMC_Flat_Interface;
using rtt_imc_test::Parser;
using rtt_imc_dd_test::IMC_DD_Interface;
using rtt_imc::Source_Builder;
using rtt_imc::Rep_Source_Builder;
using rtt_imc::DD_Source_Builder;
using rtt_imc::Opacity;
using rtt_imc::Flat_Mat_State_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Source;
using rtt_imc::Gray_Particle;
using rtt_imc::Multigroup_Particle;
using rtt_imc::Gray_Frequency;
using rtt_imc::Multigroup_Frequency;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::Comm_Patterns;
using rtt_mc::Parallel_Data_Operator;
using rtt_mc::global::soft_equiv;
using rtt_rng::Rnd_Control;
using rtt_dsxx::SP;

typedef Gray_Frequency               Gray;
typedef Gray_Particle<OS_Mesh>       GPT;
typedef SP<GPT>                      SP_GPT;

typedef Multigroup_Frequency         MG;
typedef Multigroup_Particle<OS_Mesh> MGPT;
typedef SP<MGPT>                     SP_MGPT;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// build source test for a full replication topology --> tests
// Rep_Source_Builder --> Gray_Frequency specialization

void gray_source_replication_test()
{
    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build an interface to a six cell fully replicated mesh
    SP<IMC_Flat_Interface<GPT> > interface(new IMC_Flat_Interface<GPT>(mb));

    // build a Topology: we do not use the Topology builder here because the
    // topology builder is designed to work on the host processor only -->
    // instead we will just build a Replication topology on each mesh
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));
    Parallel_Data_Operator pop(topology);

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // build a Mat_State and Opacity
    Flat_Mat_State_Builder<OS_Mesh,Gray> ob(interface);
    SP<Gray>                   frequency = ob.build_Frequency();
    SP<Mat_State<OS_Mesh> >    mat       = ob.build_Mat_State(mesh);
    SP<Opacity<OS_Mesh,Gray> > opacity   = ob.build_Opacity(mesh, frequency,
							    mat);

    if (!frequency->is_gray()) ITFAILS;

    // build a Rep_Source Builder
    Rep_Source_Builder<OS_Mesh,Gray,GPT> source_builder(
	interface, mesh, topology);

    // build the source
    SP<Source<OS_Mesh,Gray,GPT> > source = 
	source_builder.build_Source(mesh, mat, opacity, rcon, patterns);

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

    int running_rn          = 0;
    int rn_marker           = 0;

    // set the proc to get the first of the next leftover particles
    int first_leftover_proc = 0;

    // calculate census random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread  = global_ncen[cell-1] / C4::nodes();
	int leftover     = global_ncen[cell-1] - even_spread * C4::nodes();
	int shifted_node = (C4::node() + C4::nodes() - first_leftover_proc) % 
	    C4::nodes();

	rn_cen[cell-1] = running_rn + even_spread * shifted_node +
	    rtt_mc::global::min(shifted_node, leftover);

	running_rn += global_ncen[cell-1];

	// pre-combed, pre total source iteration census numbers
	local_ncen[cell-1] = even_spread;

	if (shifted_node < leftover)
	    local_ncen[cell-1]++;

	// update the proc to get the first of the next leftover particles
	first_leftover_proc = (first_leftover_proc + leftover) % C4::nodes();
    }

    // set random number id marker to first volume source id
    rn_marker = running_rn;

    // reset the proc to get the first of the next leftover particles
    first_leftover_proc = 0;

    // calculate volume random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread  = global_nvol[cell-1] / C4::nodes();
	int leftover     = global_nvol[cell-1] - even_spread * C4::nodes();
	int shifted_node = (C4::node() + C4::nodes() - first_leftover_proc) % 
	    C4::nodes();

	rn_vol[cell-1] = running_rn + even_spread * shifted_node +
	    rtt_mc::global::min(shifted_node, leftover);

	running_rn += global_nvol[cell-1];

	local_nvol[cell-1] = even_spread;

	if (shifted_node < leftover)
	    local_nvol[cell-1]++;

	// update the proc to get the first of the next leftover particles
	first_leftover_proc = (first_leftover_proc + leftover) % C4::nodes();
    }

    // add up random number ids in each cell for volume source
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	for (int np = 0; np < global_nvol[cell-1]; np++)
	    id_sum_vol[cell-1] += rn_marker + np;

	rn_marker += global_nvol[cell-1];
    }

    Check (rn_marker == running_rn);

    // reset the proc to get the first of the next leftover particles
    first_leftover_proc = 0;

    // calculate surface source random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread  = global_nss[cell-1] / C4::nodes();
	int leftover     = global_nss[cell-1] - even_spread * C4::nodes();
	int shifted_node = (C4::node() + C4::nodes() - first_leftover_proc) % 
	    C4::nodes();

	rn_ss[cell-1] = running_rn + even_spread * shifted_node +
	    rtt_mc::global::min(shifted_node, leftover);

	running_rn += global_nss[cell-1];

	local_nss[cell-1] = even_spread;

	if (shifted_node < leftover)
	    local_nss[cell-1]++;

	// update the proc to get the first of the next leftover particles
	first_leftover_proc = (first_leftover_proc + leftover) % C4::nodes();
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
	SP_GPT particle      = source->get_Source_Particle(.001);
	ss_ew               += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (cell < 1) ITFAILS;
	if (cell > 6) ITFAILS;

	if (rn_ss[cell-1] > streamid)                      ITFAILS;
	if (streamid >= rn_ss[cell-1] + local_nss[cell-1]) ITFAILS;

	// add up ids
	calc_ss_rn_sum[cell-1] += streamid;
    }

    // get volume sources
    for (int i = 0; i < source->get_nvoltot(); i++)
    {
	SP_GPT particle      = source->get_Source_Particle(.001);
	vol_ew              += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (cell < 1) ITFAILS;
	if (cell > 6) ITFAILS;

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
	SP_GPT particle      = source->get_Source_Particle(.001);
	cen_ew              += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (cell < 1) ITFAILS;
	if (cell > 6) ITFAILS;

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

    if (rtt_imc_test::passed)
	PASSMSG("Rep_Source_Builder tests pass for Gray_Frequency.");
}

//---------------------------------------------------------------------------//
// build source test for a full replication topology --> tests
// Rep_Source_Builder --> Multigroup_Frequency specialization

void mg_source_replication_test()
{
    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));
    rtt_rng::rn_stream = 0;

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build an interface to a six cell fully replicated mesh
    bool constant_group_opacities = true;
    SP<IMC_Flat_Interface<MGPT> > interface(
	new IMC_Flat_Interface<MGPT>(mb, constant_group_opacities));

    // build a Topology: we do not use the Topology builder here because the
    // topology builder is designed to work on the host processor only -->
    // instead we will just build a Replication topology on each mesh
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));
    Parallel_Data_Operator pop(topology);

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // build a Mat_State and Opacity
    Flat_Mat_State_Builder<OS_Mesh,MG> ob(interface);
    SP<MG>                   frequency = ob.build_Frequency();
    SP<Mat_State<OS_Mesh> >  mat       = ob.build_Mat_State(mesh);
    SP<Opacity<OS_Mesh,MG> > opacity   = ob.build_Opacity(mesh, frequency, mat);

    if (frequency->is_gray())        ITFAILS;
    if (!frequency->is_multigroup()) ITFAILS;

    // calculate Planck integrals
    vector<double> sigma_Planck(6, 0.0);
    vector<double> norm_Planck(6, 0.0);
    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	double T = mat->get_T(i);

	for (int g = 1; g <= frequency->get_num_groups(); g++)
	{
	    pair<double,double> bds = frequency->get_group_boundaries(g);
	    sigma_Planck[i-1]      += opacity->get_sigma_abs(i,g) * 
		rtt_cdi::CDI::integratePlanckSpectrum(bds.first, bds.second,
						      T);
	}
	
	norm_Planck[i-1] = rtt_cdi::CDI::integratePlanckSpectrum(0.01, 100.0, 
								 T);
    }

    // build a Rep_Source Builder
    Rep_Source_Builder<OS_Mesh,MG,MGPT> source_builder(
	interface, mesh, topology);

    // build the source
    SP<Source<OS_Mesh,MG,MGPT> > source = 
	source_builder.build_Source(mesh, mat, opacity, rcon, patterns);

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
    // initial random number stream ids, these are not necessarily the final
    // census numbers; however, because there is no combing, they are
    // equivalent to the final census numbers
    vector<int> global_nvol(mesh->num_cells(), 1);
    vector<int> global_ncen(mesh->num_cells(), 153);
    vector<int> global_nss(mesh->num_cells(), 0);
    {
	global_ncen[3] = 116;
	global_ncen[4] = 116;
	global_ncen[5] = 116;
	global_nvol[3] = 2;
	global_nvol[4] = 2;
	global_nvol[5] = 2;
	global_nss[0]  = 92;
	global_nss[1]  = 92;
    }

    // local source numbers per processor per species
    vector<int> local_nvol(mesh->num_cells(), 0);
    vector<int> local_ncen(mesh->num_cells(), 0);
    vector<int> local_nss(mesh->num_cells(), 0);

    // sum of ids for surface source and volume source per cell
    vector<int> id_sum_ss(mesh->num_cells(), 0);
    vector<int> id_sum_vol(mesh->num_cells(), 0);

    int running_rn          = 0;
    int rn_marker           = 0;

    // set the proc to get the first of the next leftover particles
    int first_leftover_proc = 0;

    // calculate census random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread  = global_ncen[cell-1] / C4::nodes();
	int leftover     = global_ncen[cell-1] - even_spread * C4::nodes();
	int shifted_node = (C4::node() + C4::nodes() - first_leftover_proc) % 
	    C4::nodes();

	rn_cen[cell-1] = running_rn + even_spread * shifted_node +
	    rtt_mc::global::min(shifted_node, leftover);

	running_rn += global_ncen[cell-1];

	// pre-combed, pre total source iteration census numbers
	local_ncen[cell-1] = even_spread;

	if (shifted_node < leftover)
	    local_ncen[cell-1]++;

	// update the proc to get the first of the next leftover particles
	first_leftover_proc = (first_leftover_proc + leftover) % C4::nodes();
    }

    // set random number id marker to first volume source id
    rn_marker = running_rn;

    // reset the proc to get the first of the next leftover particles
    first_leftover_proc = 0;

    // calculate volume random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread  = global_nvol[cell-1] / C4::nodes();
	int leftover     = global_nvol[cell-1] - even_spread * C4::nodes();
	int shifted_node = (C4::node() + C4::nodes() - first_leftover_proc) % 
	    C4::nodes();

	rn_vol[cell-1] = running_rn + even_spread * shifted_node +
	    rtt_mc::global::min(shifted_node, leftover);

	running_rn += global_nvol[cell-1];

	local_nvol[cell-1] = even_spread;

	if (shifted_node < leftover)
	    local_nvol[cell-1]++;

	// update the proc to get the first of the next leftover particles
	first_leftover_proc = (first_leftover_proc + leftover) % C4::nodes();
    }

    // add up random number ids in each cell for volume source
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	for (int np = 0; np < global_nvol[cell-1]; np++)
	    id_sum_vol[cell-1] += rn_marker + np;

	rn_marker += global_nvol[cell-1];
    }

    Check (rn_marker == running_rn);

    // reset the proc to get the first of the next leftover particles
    first_leftover_proc = 0;

    // calculate surface source random number id
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	int even_spread  = global_nss[cell-1] / C4::nodes();
	int leftover     = global_nss[cell-1] - even_spread * C4::nodes();
	int shifted_node = (C4::node() + C4::nodes() - first_leftover_proc) % 
	    C4::nodes();

	rn_ss[cell-1] = running_rn + even_spread * shifted_node +
	    rtt_mc::global::min(shifted_node, leftover);

	running_rn += global_nss[cell-1];

	local_nss[cell-1] = even_spread;

	if (shifted_node < leftover)
	    local_nss[cell-1]++;

	// update the proc to get the first of the next leftover particles
	first_leftover_proc = (first_leftover_proc + leftover) % C4::nodes();
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
	SP_MGPT particle     = source->get_Source_Particle(.001);
	ss_ew               += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (cell < 1) ITFAILS;
	if (cell > 6) ITFAILS;

	if (rn_ss[cell-1] > streamid)                      ITFAILS;
	if (streamid >= rn_ss[cell-1] + local_nss[cell-1]) ITFAILS;

	// add up ids
	calc_ss_rn_sum[cell-1] += streamid;
    }

    // get volume sources
    for (int i = 0; i < source->get_nvoltot(); i++)
    {
	SP_MGPT particle     = source->get_Source_Particle(.001);
	vol_ew              += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (cell < 1) ITFAILS;
	if (cell > 6) ITFAILS;

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
	SP_MGPT particle     = source->get_Source_Particle(.001);
	cen_ew              += particle->get_ew();

	int streamid = particle->get_random().get_num();
	int cell     = particle->get_cell();

	if (cell < 1) ITFAILS;
	if (cell > 6) ITFAILS;

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
    double hand_vol = 12.747095;
    double hand_ss  = 325.909753;
    double hand_cen = 1436.471652;

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
    if (fabs(evn - 0.466850) > 1.e-4 * 0.466850) ITFAILS;
    evn = source_builder.get_evol_net(2);
    if (fabs(evn - 0.466850) > 1.e-4 * 0.466850) ITFAILS;
    evn = source_builder.get_evol_net(3);
    if (fabs(evn - 0.466850) > 1.e-4 * 0.466850) ITFAILS;

    evn = source_builder.get_evol_net(4);
    if (fabs(evn - 2.620023) > 1.e-4 * 2.620023) ITFAILS;
    evn = source_builder.get_evol_net(5);
    if (fabs(evn - 2.620023) > 1.e-4 * 2.620023) ITFAILS;
    evn = source_builder.get_evol_net(6);
    if (fabs(evn - 2.620023) > 1.e-4 * 2.620023) ITFAILS;

    if (rtt_imc_test::passed)
	PASSMSG("Rep_Source_Builder tests pass for Multigroup_Frequency.");
}

//---------------------------------------------------------------------------//
// test the source builder for a fully domain-decomposed (DD) mesh. -->
// Gray_Frequency specialization

void gray_source_DD_test()
{
    // the DD test mesh is only for four processors
    if (C4::nodes() != 4)
	return;

    // build a random number controller and reinitialize the global random
    // number stream number, which was incremented in previous source_builder
    // tests.
    SP<Rnd_Control> rcon(new Rnd_Control(347223));
    rtt_rng::rn_stream = 0;

    // build each processor's portion of the DD Mesh 
    // (9-cell DD mesh from IMC_DD_Test.hh)
    SP<OS_Mesh> mesh = rtt_imc_dd_test::build_DD_Mesh(); 

    // check number of cells in mesh on each processor
    if (!((C4::node() != 3) ? (mesh->num_cells() == 2) : true)) ITFAILS;
    if (!((C4::node() == 3) ? (mesh->num_cells() == 3) : true)) ITFAILS;

    // build the topology for the DD test mesh
    SP<Topology> topology = rtt_imc_dd_test::build_DD_Topology();
    Parallel_Data_Operator pop(topology);

    // build a comm_patterns
    SP<Comm_Patterns> patterns(new Comm_Patterns());
    patterns->calc_patterns(topology);

    // check parallel scheme
    if (topology->get_parallel_scheme() != "DD") ITFAILS;

    // check topology's value of num_cells() compared to mesh's value
    if (mesh->num_cells() != topology->num_cells(C4::node())) ITFAILS;

    // need an interface object (not explicitly used in
    // calc_fullDD_rn_fields).
    SP<IMC_DD_Interface<GPT> > interface(new IMC_DD_Interface<GPT>(
					     mesh->num_cells())); 

    // check one thing from the interface
    if (interface->get_delta_t() != 0.001) ITFAILS;

    // build a Mat_State and Opacity
    Flat_Mat_State_Builder<OS_Mesh,Gray> ob(interface);
    SP<Gray>                   frequency = ob.build_Frequency();
    SP<Mat_State<OS_Mesh> >    mat       = ob.build_Mat_State(mesh);
    SP<Opacity<OS_Mesh,Gray> > opacity   = ob.build_Opacity(mesh, frequency, 
							    mat);

    // build a DD_Source Builder
    DD_Source_Builder<OS_Mesh,Gray,GPT> source_builder(
	interface, mesh, topology);

    // build the source
    SP<Source<OS_Mesh,Gray,GPT> > source = 
	source_builder.build_Source(mesh, mat, opacity, rcon, patterns); 

    // <<<<<< SET REFERENCE VARIABLES >>>>>>

    // set reference variables
    vector<double>        temp(mesh->num_cells(), 0.0);
    vector<double>        dens(mesh->num_cells(), 0.0); 
    vector<double>        kapp(mesh->num_cells(), 0.0); 
    vector<double>        kapt(mesh->num_cells(), 0.0);
    vector<double>         c_v(mesh->num_cells(), 0.0); 
    vector<double>     evolext(mesh->num_cells(), 100.0);
    vector<double>     rad_src(mesh->num_cells(), 200.0);
    vector<double>       rtemp(mesh->num_cells(), 10.0);
    vector<double>    ref_evol(mesh->num_cells(), 0.0);
    vector<double> ref_evolnet(mesh->num_cells(), 0.0);
    vector<double>    ref_ecen(mesh->num_cells(), 0.0);
    vector<double> ref_matvsrc(mesh->num_cells(), 0.0);
    double radconst       = 0.01372;
    double sol            = 299.792458;
    double ref_evoltot    = 0.0;
    double ref_esstot     = sol*radconst/4.0 *20*20*20*20 * 0.001;
    double ref_ecentot    = 0.0;
    double ref_matvsrctot = 0.0;

    for (int c = 0; c < mesh->num_cells(); c++)
    {
	temp[c] = 3.0*(C4::node()+c+1);
	dens[c] = C4::node() + c + 1;
	kapp[c] = 2.0*C4::node() + c + 1;
	kapt[c] = 2.0*(C4::node()+c+1);
	c_v[c]  = 3.0*C4::node() + c + 1;
	double fleck    = 1.0/(1.0 +  4.0 * radconst * temp[c]*temp[c]*temp[c] * 
			       sol * 0.001 * kapp[c] / c_v[c]);
	ref_evolnet[c]  = fleck * kapp[c] * dens[c]* radconst * sol *
	                   (temp[c]*temp[c]*temp[c]*temp[c]) * 0.001; 
	// ext src
	ref_evol[c]     =  ref_evolnet[c] + 100.0 * (1.0-fleck) * 0.001;
	// add rad_source to evol
	ref_evol[c]    += 200.0 * 0.001; 
	ref_matvsrc[c]  = 100.0 * fleck * 0.001;
	ref_matvsrctot += ref_matvsrc[c];
	ref_ecen[c]     = radconst * 1.0e4;
	ref_evoltot    += ref_evol[c];
	ref_ecentot    += ref_ecen[c];
    }
	
    // <<<<<< CHECK LOCAL VOLUME EMISSION ENERGY >>>>>>>> 

    // check net volume emission (each proc has 2 cells, except proc 3 has 3)
    if (mesh->num_cells() > 0)
	if (!soft_equiv(ref_evolnet[0],source_builder.get_evol_net(1)))
	    ITFAILS;
    if (mesh->num_cells() > 1)
	if (!soft_equiv(ref_evolnet[1],source_builder.get_evol_net(2)))
	    ITFAILS; 
    if (mesh->num_cells() > 2)
	if (!soft_equiv(ref_evolnet[2],source_builder.get_evol_net(3)))
	    ITFAILS; 

    // do processor-dependent evol checks only if there is no global energy
    // loss
    double evol_loss     = source_builder.get_eloss_vol();
    double check_evoltot = source_builder.get_evoltot();
    if (soft_equiv(evol_loss, 0.0, 1.0e-12))
	if (!soft_equiv(check_evoltot, ref_evoltot, 1.0e-12)) ITFAILS;

    // do global evol check
    C4::gsum(ref_evoltot);	
    C4::gsum(check_evoltot);
    if (!soft_equiv(check_evoltot+evol_loss, ref_evoltot, 1.0e-12)) ITFAILS;


    // <<<<<< SURFACE SOURCE ENERGY >>>>>>>>

    // do processor-dependent ess checks only if there is no global energy loss
    double ess_loss = source_builder.get_eloss_ss();
    if (soft_equiv(ess_loss, 0.0, 1.0e-12))
    {
	// each processor has one surface source--compare to hand calc
	if (!soft_equiv(source_builder.get_esstot(), 164.526101, 1.0e-6))
	    ITFAILS; 
	// each processor has one surface source--compare to reference calc
	if (!soft_equiv(source_builder.get_esstot(), ref_esstot)) ITFAILS; 
    }

    // do global ess checks (energies same on each of 4 processors)
    double global_esstot = source_builder.get_esstot();
    C4::gsum(global_esstot);

    if (!soft_equiv(global_esstot+ess_loss, 4.0*ref_esstot, 1.0e-12))
	ITFAILS;  


    // <<<<<< MATERIAL VOLUME SOURCE ENERGY >>>>>>>>

    // check material volume source
    if (mesh->num_cells() > 0)
	if (!soft_equiv(ref_matvsrc[0],source_builder.get_mat_vol_src(1)))
	    ITFAILS;
    if (mesh->num_cells() > 1)
	if (!soft_equiv(ref_matvsrc[1],source_builder.get_mat_vol_src(2)))
	    ITFAILS; 
    if (mesh->num_cells() > 2)
	if (!soft_equiv(ref_matvsrc[2],source_builder.get_mat_vol_src(3)))
	    ITFAILS; 

    // check on-processor totals
    if (!soft_equiv(ref_matvsrctot, source_builder.get_mat_vol_srctot()))
	ITFAILS; 


    // <<<<<< INITIAL CENSUS ENERGY >>>>>>>>

    // global census energy check
    double global_ecentot = source_builder.get_initial_census_energy();
    double ref_global_ecentot = ref_ecentot;
    C4::gsum(ref_global_ecentot);
    double check_ecentot = global_ecentot + source_builder.get_eloss_cen(); 
    if (!soft_equiv(check_ecentot, ref_global_ecentot, 1.0e-12))     ITFAILS;


    // <<<<<< NUMBER OF VOLUME EMISSION PARTICLES >>>>>>

    // check consistency between local cell values and global total
    vector<int> nvol(mesh->num_cells(),0);
    int nvoltot = 0;
    for (int c = 1; c <= mesh->num_cells(); c++)
    {
	nvol[c-1] = source_builder.get_nvol(c);
	nvoltot += nvol[c-1];
    }
    int global_nvoltot = nvoltot;
    C4::gsum(global_nvoltot);
    if (global_nvoltot != source_builder.get_nvoltot()) ITFAILS;


    // <<<<<< GLOBAL TOTAL NUMBER OF PARTICLES >>>>>>>> 

    // contrived (not hand checked) check on global total number of particles
    int global_nsstot  = source_builder.get_nsstot();
    int global_ncentot = source_builder.get_ncentot();
    if (global_nvoltot + global_nsstot + global_ncentot != 991) ITFAILS; 


    // <<<<<< CALC VOLUME EMISSION RANDOM NUMBER STREAM INFO >>>>>>

    // global vector of number of volume emission particles
    int *global_nvol = new int[9];
    for (int i = 0; i < 9; i++)
	global_nvol[i] = 0;

    // loop over local cells and map local nvol to global vector
    for (int lc = 1; lc <= mesh->num_cells(); lc++)
    {
	int gc = 2*C4::node() + lc;
	global_nvol[gc-1] = nvol[lc-1];
    }
    C4::gsum(global_nvol,9);
    
    // calc volrn and construct nvol rn checksum
    vector<int> volrn(9, 0);
    int checksum_vol_rn = 0;
    int offset = 0;
    for (int gc = 0; gc < 9; gc++)
    {
	volrn[gc] = global_ncentot + offset;

	for (int np = 0; np < global_nvol[gc]; np++)
	    checksum_vol_rn += volrn[gc] + np;

	offset += global_nvol[gc];
    }

    // <<<<<< RETRIEVE PARTICLES FROM SOURCE >>>>>>

    // sums of surface source and vol emiss ew's
    double calc_ss_ewtot  = 0.0;    
    double calc_vol_ewtot = 0.0;
    double calc_cen_ewtot = 0.0;

    // sum of volume emission particle random number stream id's
    int sum_vol_rn = 0;

    // get surface sources; weakly check random number stream id's
    for (int i = 0; i < source->get_nsstot(); i++)
    {
	SP_GPT particle      = source->get_Source_Particle(.001);
	calc_ss_ewtot       += particle->get_ew();

	int ssrnid = particle->get_random().get_num();
	if (ssrnid < global_ncentot + global_nvoltot)                ITFAILS; 
	if (ssrnid > 990)                                            ITFAILS; 
    }
    C4::gsum(calc_ss_ewtot);

    if (!soft_equiv(calc_ss_ewtot + ess_loss, 4.0*ref_esstot, 1.0e-12))
	ITFAILS; 

    // get volume sources; check random number stream id's
    for (int i = 0; i < source->get_nvoltot(); i++)
    {
	SP_GPT particle = source->get_Source_Particle(.001);
	calc_vol_ewtot += particle->get_ew();

	int volrnid     = particle->get_random().get_num();
	sum_vol_rn     += volrnid;
	int cell        = particle->get_cell();
	int global_cell = 2*C4::node() + cell;

	int lo_id_limit = volrn[global_cell-1];
	int hi_id_limit = volrn[global_cell-1]+global_nvol[global_cell-1];

	if (volrnid <  lo_id_limit)                                  ITFAILS;
	if (volrnid >= hi_id_limit)                                  ITFAILS;
    }
    C4::gsum(calc_vol_ewtot);
    C4::gsum(sum_vol_rn);

    if (!soft_equiv(calc_vol_ewtot+evol_loss, ref_evoltot, 1.0e-12)) ITFAILS;
    if (sum_vol_rn != checksum_vol_rn)                               ITFAILS; 


    // get census sources; weakly check random number stream id's
    for (int i = 0; i < source->get_ncentot(); i++)
    {
	SP_GPT particle      = source->get_Source_Particle(.001);
	calc_cen_ewtot      += particle->get_ew();
	int lcell            = particle->get_cell();
	int gcell            = topology->global_cell(lcell);
	if (gcell < 2*C4::node() + 1)                                ITFAILS;
	if (C4::node() < 3)
	    if (gcell > 2*(C4::node()+1))                            ITFAILS;
	if (C4::node() == 3)
	    if (gcell > 2*(C4::node()+1) + 1)                        ITFAILS;
	
	int cenrnid  = particle->get_random().get_num();
	if (cenrnid >= global_ncentot)                               ITFAILS;
	if (cenrnid < 0)                                             ITFAILS;
    }
    C4::gsum(calc_cen_ewtot);

    if (!soft_equiv(calc_cen_ewtot + source_builder.get_eloss_cen(),
		    ref_global_ecentot, 1.0e-9))                     ITFAILS; 

    // the source should be empty
    if (*source) ITFAILS;

    if (rtt_imc_test::passed)
	PASSMSG("DD_Source_Builder tests pass for Gray_Frequency.");
}

//---------------------------------------------------------------------------//

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

    try
    {
	// >>> UNIT TESTS

	// full replication source test for gray and multigroup
	gray_source_replication_test();
	mg_source_replication_test();

	// source builder test on full domain decomposition; we only need to
	// test gray because all of the Frequency dependencies are in the
	// Source_Builder base class
	gray_source_DD_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstSource_Builder, " << ass.what()
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
	    cout << "**** tstSource_Builder Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstSource_Builder on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstSource_Builder.cc
//---------------------------------------------------------------------------//
