//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstMesh_Operations.cc
 * \author Thomas M. Evans
 * \date   Mon Dec 20 16:44:16 1999
 * \brief  Mesh_Operations test class
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Mesh_Operations.hh"
#include "../Opacity_Builder.hh"
#include "../Opacity.hh"
#include "../Mat_State.hh"
#include "../Release.hh"
#include "mc/Topology.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

using namespace std;

using rtt_imc_test::IMC_Interface;
using rtt_imc::Mesh_Operations;
using rtt_imc::Opacity_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Opacity;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_dsxx::SP;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

//---------------------------------------------------------------------------//
// Test T4 slope tilt

void T4_slope_test()
{
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

    // build a Mesh_Operations class for OS_Mesh
    Mesh_Operations<OS_Mesh> mesh_op(mesh, mat, topology);

    // get the t4 slopes
    OS_Mesh::CCVF<double> slopes = mesh_op.get_t4_slope();

    // test the t4 slope values
    if (slopes.size() != 2)                  ITFAILS;
    if (slopes.size(1) != mesh->num_cells()) ITFAILS;
    if (slopes.size(2) != mesh->num_cells()) ITFAILS;

    if (slopes(1,1) != 0) ITFAILS;
    if (slopes(1,2) != 0) ITFAILS;
    if (slopes(1,3) != 0) ITFAILS;
    if (slopes(1,4) != 0) ITFAILS;
    if (slopes(1,5) != 0) ITFAILS;
    if (slopes(1,6) != 0) ITFAILS;

    if (slopes(2,1) != 10000) ITFAILS;
    if (slopes(2,2) != 10000) ITFAILS;
    if (slopes(2,3) != 10000) ITFAILS;
    if (slopes(2,4) != 75000) ITFAILS;
    if (slopes(2,5) != 75000) ITFAILS;
    if (slopes(2,6) != 75000) ITFAILS;

    // test sampling

    // make two identical random number objects
    Rnd_Control rcon(39567);
    Sprng ran1 = rcon.get_rn(10);
    Sprng ran2 = rcon.get_rn(10);
    if (ran1.ran() != ran2.ran()) ITFAILS;

    // sample coordinate position in cell 5
    vector<double> reference_r(2);
    double T4 = 20 * 20 * 20 * 20;
   
    // x coordinate position
    double bx    = T4 - slopes(1,5) * .5;
    double probx = .5 * bx / T4;

    if (ran1.ran() <= probx)
	reference_r[0] = 1.0 - (1.0) * sqrt(ran1.ran());
    else
	reference_r[0] = 0.0 + (1.0) * sqrt(ran1.ran());

    // y coordinate position
    double by    = T4 - slopes(2,5) * (2.0) * .5;
    double proby = .5 * by / T4;

    if (ran1.ran() <= proby)
	reference_r[1] = 3.0 - (2.0) * sqrt(ran1.ran());
    else
	reference_r[1] = 1.0 + (2.0) * sqrt(ran1.ran());

    // get sampling from mesh operations
    vector<double> r = mesh_op.sample_pos_tilt(5, mat->get_T(5), ran2);

    // check for equivalence
    if (reference_r != r) ITFAILS;
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

    // run test on each processor
    T4_slope_test();
	

    // status of test
    cout << endl;
    cout <<     "*******************************************" << endl; 
    if (passed) 
    {
        cout << "**** Mesh_Operations Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "*******************************************" << endl;
    cout << endl;

    cout << "Done testing Mesh_Operations on node: " << C4::node() << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstMesh_Operations.cc
//---------------------------------------------------------------------------//
