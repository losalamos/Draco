//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstSource_Builder.cc
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
#include "mc/Rep_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

using rtt_imc_test::IMC_Interface;
using rtt_imc::Source_Builder;
using rtt_imc::Rep_Source_Builder;
using rtt_imc::Opacity;
using rtt_imc::Opacity_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Source;
using rtt_imc::Topology_Builder;
using rtt_imc::Particle;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_rng::Rnd_Control;
using dsxx::SP;

typedef SP<Particle<OS_Mesh> > SP_Particle;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

//---------------------------------------------------------------------------//

void source_replication_test()
{
    // build a random number controller
    SP<Rnd_Control> rcon(new Rnd_Control(347223));

    // build an interface to a six cell fully replicated mesh
    SP<IMC_Interface> interface(new IMC_Interface);

    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    OS_Builder mb(interface);
    SP<OS_Mesh> mesh = mb.build_Mesh();

    // build a Topology: we do not use the Topology builder here because the
    // topology builder is designed to work on the host processor only -->
    // instead we will just build a Replication topology on each mesh
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a Mat_State and Opacity
    Opacity_Builder<OS_Mesh> ob(interface);
    SP<Mat_State<OS_Mesh> > mat    = ob.build_Mat(mesh);
    SP<Opacity<OS_Mesh> > opacity  = ob.build_Opacity(mesh, mat);

    // build a Rep_Source Builder
    Rep_Source_Builder<OS_Mesh> source_builder(interface, mesh, topology);

    // build the source
    SP<Source<OS_Mesh> > source = source_builder.build_Source(mesh, mat,
							      opacity, rcon);

    // get particles out of the source and add up the energy weights
    double cen_ew = 0.0;
    double vol_ew = 0.0;
    double ss_ew  = 0.0;
    
    // get surface sources
    for (int i = 0; i < source->get_nsstot(); i++)
    {
	SP_Particle particle = source->get_Source_Particle(.001);
	ss_ew               += particle->get_ew();
    }

    // get volume sources
    for (int i = 0; i < source->get_nvoltot(); i++)
    {
	SP_Particle particle = source->get_Source_Particle(.001);
	vol_ew              += particle->get_ew();
    }

    // get census sources
    for (int i = 0; i < source->get_ncentot(); i++)
    {
	SP_Particle particle = source->get_Source_Particle(.001);
	cen_ew              += particle->get_ew();
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
