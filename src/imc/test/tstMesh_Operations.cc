//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstMesh_Operations.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:04:34 2001
 * \brief  Mesh_Operations test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "DD_Mesh.hh"
#include "../Mesh_Operations.hh"
#include "../Opacity_Builder.hh"
#include "../Opacity.hh"
#include "../Mat_State.hh"
#include "../Release.hh"
#include "../Global.hh"
#include "mc/Topology.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/Comm_Patterns.hh"
#include "mc/RZWedge_Mesh.hh"
#include "mc/Sampler.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

using namespace std;

using rtt_imc_test::IMC_Interface;
using rtt_imc_test::Parser;
using rtt_imc::Mesh_Operations;
using rtt_imc::Opacity_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Opacity;
using rtt_mc::global::soft_equiv;
using rtt_mc::sampler::sample_general_linear;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::Comm_Patterns;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// Test T4 slope tilt

void T4_slope_test()
{
    // build a FULL mesh --> this mesh will be fully replicated on all
    // processors in the test
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // build an interface to a six cell fully replicated mesh
    SP<IMC_Interface> interface(new IMC_Interface(mb));

    // build a Topology: we do not use the Topology builder here because the
    // topology builder is designed to work on the host processor only -->
    // instead we will just build a Replication topology on each mesh
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // build a Mat_State and Opacity
    Opacity_Builder<OS_Mesh> ob(interface);
    SP<Mat_State<OS_Mesh> > mat    = ob.build_Mat(mesh);
    SP<Opacity<OS_Mesh> > opacity  = ob.build_Opacity(mesh, mat);

    // Build a Comm_Patterns
    SP<Comm_Patterns> comm_patterns(new Comm_Patterns());
    comm_patterns->calc_patterns(topology);

    // build a Mesh_Operations class for OS_Mesh
    Mesh_Operations<OS_Mesh> mesh_op(mesh, mat, topology, comm_patterns);

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

void T4_slope_test_DD()
{
    if (C4::nodes() != 4)
	return;

    // local meshes 9 cell (4 cells on processor 0, 5 cells on processor 1)
    SP<OS_Mesh>             mesh;
    SP<Topology>            topology;
    SP<Mat_State<OS_Mesh> > mat;

    mesh     = rtt_imc_test::build_Mesh();
    topology = rtt_imc_test::build_Topology();
    mat      = rtt_imc_test::build_Mat(mesh);

    // make and calculate comm_patterns
    SP<Comm_Patterns> comm_patterns(new Comm_Patterns());
    comm_patterns->calc_patterns(topology);

    // build a Mesh_Operations class for OS_Mesh
    Mesh_Operations<OS_Mesh> mesh_op(mesh, mat, topology, comm_patterns);

    // get the slopes
    OS_Mesh::CCVF<double> slopes = mesh_op.get_t4_slope();
    if (slopes.size() != 2) ITFAILS;

    // now check out the data
    if (C4::node() == 0)
    {
	if (slopes.size(1) != 2)   ITFAILS;
	if (slopes(1,1) != 0.0)    ITFAILS;
	if (slopes(1,2) != 2667.5) ITFAILS;
	if (slopes.size(2) != 2)   ITFAILS;
	if (slopes(2,1) != 5335)   ITFAILS;
	if (slopes(2,2) != 12259)  ITFAILS;
    }
    else if (C4::node() == 1)
    {
	if (slopes.size(1) != 2)   ITFAILS;
	if (slopes(1,1) != 5335)   ITFAILS;
	if (slopes(1,2) != 6924)   ITFAILS;
	if (slopes.size(2) != 2)   ITFAILS;
	if (slopes(2,1) != 6924)   ITFAILS;
	if (slopes(2,2) != 10530)  ITFAILS;
    }
    else if (C4::node() == 2)
    {
	if (slopes.size(1) != 2)   ITFAILS;
	if (slopes(1,1) != 3462)   ITFAILS;
	if (slopes(1,2) != 0)      ITFAILS;
	if (slopes.size(2) != 2)   ITFAILS;
	if (slopes(2,1) != 10530)  ITFAILS;
	if (slopes(2,2) != 7862.5) ITFAILS;
    }
    else if (C4::node() == 3)
    {
	if (slopes.size(1) != 3)   ITFAILS;
	if (slopes(1,1) != 0)      ITFAILS;
	if (slopes(1,2) != 0)      ITFAILS;
	if (slopes(1,3) != 0)      ITFAILS;
	if (slopes.size(2) != 3)   ITFAILS;
	if (slopes(2,1) != 0)      ITFAILS;
	if (slopes(2,2) != 0)      ITFAILS;
	if (slopes(2,3) != 0)      ITFAILS;
    }

    // we already checked sample_pos_tilt above, since the slopes are correct 
    // we can assume that sample_pos_tilt is also correct
}

//===========================================================================//
// RZWEDGE_MESH SPECIALIZATION
//===========================================================================//

void T4_slope_test_AMR()
{
    // make an AMR mesh
    SP<RZWedge_Mesh> mesh = rtt_imc_test::make_RZWedge_Mesh_AMR(30.0);

    // make a full replication topology
    SP<Topology> topology(new Rep_Topology(mesh->num_cells()));

    // make a Mat_State
    RZWedge_Mesh::CCSF<double> temps(mesh);
    for (int i = 1; i <= 3; i++)
	temps(i) = 1.0;
    for (int i = 4; i <= 6; i++)
	temps(i) = 2.0;
    for (int i = 7; i <= 9; i++)
	temps(i) = 3.0;
    SP<Mat_State<RZWedge_Mesh> > mat(new Mat_State<RZWedge_Mesh>
				     (temps, temps, temps, temps));

    // make a comm_patterns
    SP<Comm_Patterns> cp(new Comm_Patterns());
    cp->calc_patterns(topology);
    
    // Make Mesh_Operations instance
    Mesh_Operations<RZWedge_Mesh> mesh_op(mesh, mat, topology, cp);

    // test values of slopes
    RZWedge_Mesh::CCVF<double> slopes = mesh_op.get_t4_slope();
    if (slopes.size() != 3)                  ITFAILS;
    if (slopes.size(1) != mesh->num_cells()) ITFAILS;
    if (slopes.size(2) != mesh->num_cells()) ITFAILS;
    if (slopes.size(3) != mesh->num_cells()) ITFAILS;

    for (int i = 1; i <= mesh->num_cells(); i++)
	if (!soft_equiv(slopes(2, i), 0.)) ITFAILS;

    // check values
    if (!soft_equiv(slopes(1, 1), 0.))             ITFAILS;
    if (!soft_equiv(slopes(1, 2), 0.))             ITFAILS;
    if (!soft_equiv(slopes(1, 3), 0.))             ITFAILS;
    if (!soft_equiv(slopes(1, 4), 15.))            ITFAILS;
    if (!soft_equiv(slopes(1, 5), 15.))            ITFAILS;
    if (!soft_equiv(slopes(1, 6), 64.))            ITFAILS;
    if (!soft_equiv(slopes(1, 7), 53.3333, 1e-4))  ITFAILS;
    if (!soft_equiv(slopes(1, 8), 130.0))          ITFAILS;
    if (!soft_equiv(slopes(1, 9), 0.))             ITFAILS;

    if (!soft_equiv(slopes(3, 1), 0.))             ITFAILS;
    if (!soft_equiv(slopes(3, 2), 0.))             ITFAILS;
    if (!soft_equiv(slopes(3, 3), 0.))             ITFAILS;
    if (!soft_equiv(slopes(3, 4), 0.))             ITFAILS;
    if (!soft_equiv(slopes(3, 5), 21.6667, 1e-4))  ITFAILS;
    if (!soft_equiv(slopes(3, 6), 64.))            ITFAILS;
    if (!soft_equiv(slopes(3, 7), 130.0))          ITFAILS;
    if (!soft_equiv(slopes(3, 8), 43.3333, 1e-4))  ITFAILS;
    if (!soft_equiv(slopes(3, 9), 0.))             ITFAILS;

    // test sampling in cell 6

    // make two identical random number objects
    Rnd_Control rcon(39567);
    Sprng ran1 = rcon.get_rn(10);
    Sprng ran2 = rcon.get_rn(10);
    if (ran1.ran() != ran2.ran()) ITFAILS;

    // sample by hand
    vector<double> ref(3);
    {    
	// get x-dimension cell extents
	double half_delx = 0.5 * (1.5 - 1.0);

	// calculate weighting function values at x-dimension cell extents
	double half_delw = 64 * half_delx;
	double low_w     = 16 - half_delw;
	double high_w    = 16 + half_delw;

	// calculate posititve y-values corresponding to x-dimension cell
	// extents
	double loy = 1.0 * .267949;
	double hiy = 1.5 * .267949;

	// calc extents of function (wt-fn * half_y_value) from which to
	// sample
	double lof = loy * low_w;
	double hif = hiy * high_w;
	
	// sample x
	ref[0] = sample_general_linear(ran1, 1.0, 1.5, lof, hif);

	// calculate y-values corresponding to sampled x-value
	double pos_y = ref[0] * .267949;
	double neg_y = -pos_y;

	// sample y uniformly
	ref[1] = neg_y + ran1.ran() * (pos_y - neg_y);

	// get z-dimension cell extents
	double half_delz = 0.5 * (2.5 - 2.0);
	
	// calculate weighting function at z-dimension cell extents
	half_delw = 64 * half_delz;
	low_w     = 16 - half_delw;
	high_w    = 16 + half_delw;

	// sample z
	ref[2] = sample_general_linear(ran1, 2.0, 2.5, low_w, high_w);
    }
    
    // sample position using mesh_op
    vector<double> r = mesh_op.sample_pos_tilt(6, mat->get_T(6), ran2);
    if (r.size() != ref.size()) ITFAILS;
    
    for (int i = 0; i < 2; i++)
	if (!soft_equiv(r[i], ref[i], 1e-6)) ITFAILS;
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

	// run test on each processor for OS_Mesh
	T4_slope_test();
	T4_slope_test_DD();

	// run test on each processor for RZWedge_Mesh
	T4_slope_test_AMR();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstMesh_Operations, " << ass.what()
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
	    cout << "**** tstMesh_Operations Test: PASSED on" 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstMesh_Operations on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstMesh_Operations.cc
//---------------------------------------------------------------------------//
