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
#include "mc/Rep_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"

#include <vector>
#include <string>
#include <iostream>

using namespace std;

using rtt_imc_test::IMC_Interface;
using rtt_imc::Source_Builder;
using rtt_imc::Rep_Source_Builder;
using rtt_imc::Opacity;
using rtt_imc::Opacity_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Source;
using rtt_imc::Topology_Builder;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_rng::Rnd_Control;
using dsxx::SP;

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

    cout << "evoltot        " << source_builder.get_evoltot() << endl;
    cout << "mat_vol_srctot " << source_builder.get_mat_vol_srctot() << endl;
    cout << "esstot         " << source_builder.get_esstot() << endl;
    cout << "eloss_vol      " << source_builder.get_eloss_vol() << endl;
    cout << "eloss_ss       " << source_builder.get_eloss_ss() << endl;
    cout << "eloss_cen      " << source_builder.get_eloss_cen() << endl;

    cout << "nvoltot        " << source_builder.get_nvoltot() << endl;
    cout << "nsstot         " << source_builder.get_nsstot() << endl;
    cout << "ncentot        " << source_builder.get_ncentot() << endl;

    cout << *source << endl;
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
