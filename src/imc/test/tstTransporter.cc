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
#include "../Particle_Buffer.hh"
#include "../Particle.hh"
#include "../Release.hh"
#include "mc/Rep_Topology.hh"
#include "mc/OS_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Parallel_Data_Operator.hh"
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
using rtt_imc::Particle_Buffer;
using rtt_imc::Particle;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::Parallel_Data_Operator;
using rtt_rng::Rnd_Control;
using rtt_dsxx::SP;

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
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	cout << "Test: assertion failure at line " 
	     << ass.what() << endl;
	return 1;
    }
    catch(...)
    {
	cout << "HELP ME" << endl;
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
