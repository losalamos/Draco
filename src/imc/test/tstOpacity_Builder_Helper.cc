//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstOpacity_Builder_Helper.cc
 * \author Thomas M. Evans
 * \date   Mon Aug 11 18:03:47 2003
 * \brief  Opacity_Builder_Helper test.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Opacity_Builder_Helper.hh"
#include "../Opacity.hh"
#include "../Diffusion_Opacity.hh"
#include "../Release.hh"
#include "../Frequency.hh"
#include "../Mat_State.hh"
#include "../Fleck_Factors.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/Constants.hh"
#include "cdi/CDI.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <utility>

using namespace std;

using rtt_imc_test::Parser;
using rtt_imc::Opacity_Builder_Helper;
using rtt_imc::Fleck_Factors;
using rtt_imc::Mat_State;
using rtt_imc::Opacity;
using rtt_imc::Diffusion_Opacity;
using rtt_mc::OS_Builder;
using rtt_mc::OS_Mesh;
using rtt_cdi::CDI;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::OS_Mesh               MESH;
typedef rtt_imc::Gray_Frequency       G;
typedef rtt_imc::Multigroup_Frequency MG;

typedef MESH::CCSF<double>          ccsf;
typedef MESH::CCSF<vector<double> > ccvf;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

SP<Mat_State<MESH> > make_mat(SP<MESH> mesh)
{
    SP<Mat_State<MESH> > mat;
    
    ccsf density(mesh);
    ccsf temp(mesh);
    ccsf cv(mesh);

    for (int c = 1; c <= density.size(); c++)
    {
	density(c) = c + 0.1;
	temp(c)    = 2.0;
	cv(c)      = 0.1;
    }

    mat = new Mat_State<MESH>(density, temp, cv);
    return mat;
}

//---------------------------------------------------------------------------//

void gray_test()
{
    // build a parser and mesh
    SP<Parser>     parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<MESH>       mesh = mb->build_Mesh();
    if (mesh->num_cells() != 6) ITFAILS;

    // make mat state
    SP<Mat_State<MESH> > mat = make_mat(mesh);
    
    // make gray absorption opacities
    ccsf abs(mesh);
    for (int i = 1; i <= 6; i++)
	abs(i) = 212.0 + i;

    double dt  = 0.01;
    double imp = 1.0;

    // now make fleck factors
    SP<Fleck_Factors<MESH> > fleck = 
	Opacity_Builder_Helper<MESH,G>::build_Fleck_Factors(
	    mesh, mat, abs, dt, imp);

    // check fleck factors 
    if (!fleck)                   ITFAILS; 
    if (fleck->fleck.size() != 6) ITFAILS;
    
    double beta = 0.0;
    double f    = 0.0;
    for (int c = 1; c <= 6; c++)
    {
	beta = 4.0 * rtt_mc::global::a * 8.0 / (0.1 * mat->get_rho(c));
	f    = 1.0 / (1.0 + imp * beta * rtt_mc::global::c * dt * abs(c));

	if (!soft_equiv(fleck->fleck(c), f)) ITFAILS;
    }
    
    if (rtt_imc_test::passed)
    {
	ostringstream m;
	m << "Opacity_Builder_Helper ok for Gray_Frequency specialization.";
	PASSMSG(m.str());
    }
}

//---------------------------------------------------------------------------//

void mg_test()
{
    // build a parser and mesh
    SP<Parser>     parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<MESH>       mesh = mb->build_Mesh();
    if (mesh->num_cells() != 6) ITFAILS;

    // make mat state
    SP<Mat_State<MESH> > mat = make_mat(mesh);

    // make frequency
    vector<double> grps(3);
    grps[0] = 1.e-6;
    grps[1] = 0.01; 
    grps[2] = 100.0;
    
    SP<MG> frequency(new MG(grps));

    double dt  = 0.01;
    double imp = 1.0;
    
    // make opacities
    ccvf abs(mesh);
    ccvf sct(mesh);
    {
	for (int i = 1; i <= 6; i++)
	{
	    abs(i).resize(2);
	    sct(i).resize(2);
	    for (int g = 0; g < 2; g++)
	    {
		abs(i)[g] = g + 1.0;
		sct(i)[g] = 0.01;
	    }
	}
    }
    
    // build opacity and diff opacity object
    pair<SP<Opacity<MESH,MG> >, SP<Diffusion_Opacity<MESH> > > opacities;
    opacities = Opacity_Builder_Helper<MESH,MG>::build_Opacity(
	mesh, frequency, mat, abs, sct, dt, imp, true);

    if (!opacities.first)  ITFAILS;
    if (!opacities.second) ITFAILS;

    // check opacities
    {
	SP<Opacity<MESH,MG> > o = opacities.first;
	if (o->get_Frequency()->get_num_groups() != 2) ITFAILS;

	SP<Diffusion_Opacity<MESH> > d = opacities.second;
	
	for (int i = 1; i <= 6; i++)
	{
	    double b = 0.0;
	    double r = 0.0;
	    pair<double,double> bounds;

	    double bsum  = 0.0;
	    double rsum  = 0.0;
	    double bssum = 0.0;
	    double rssum = 0.0;
	    double beta  = 0.0;
	    double f     = 0.0;

	    vector<double> egc = o->get_emission_group_cdf(i);
	    if (egc.size() != 2) ITFAILS;

	    for (int g = 1; g <= 2; g++)
	    {
		if (o->get_sigma_abs(i,g) != g)        ITFAILS;
		if (o->get_sigma_thomson(i,g) != 0.01) ITFAILS;

		bounds = o->get_Frequency()->get_group_boundaries(g);
		CDI::integrate_Rosseland_Planckian_Spectrum(bounds.first,
							    bounds.second,
							    2.0, b, r);

		rsum += r;
		bsum += b;

		rssum += r * 1.0 / (g);
		bssum += b * g;

		if (!soft_equiv(egc[g-1], bssum)) ITFAILS;
	    }

	    if (!soft_equiv(o->get_integrated_norm_Planck(i), bsum))
		ITFAILS;
	    if (!soft_equiv(o->get_integrated_sigma_times_Planck(i), bssum))
		ITFAILS;

	    double planck = o->get_integrated_sigma_times_Planck(i) / 
		o->get_integrated_norm_Planck(i);
	    
	    beta = 4.0 * rtt_mc::global::a * 8.0 / (0.1 * mat->get_rho(i));
	    f    = 1.0 / (1.0 + imp * beta * rtt_mc::global::c * dt * planck);
	    
	    if (!soft_equiv(o->get_fleck(i), f))                      ITFAILS;
	    if (!soft_equiv(d->get_Rosseland_opacity(i), rsum/rssum)) ITFAILS;
	}
    }

    if (rtt_imc_test::passed)
    {
	ostringstream m;
	m << "Opacity_Builder_Helper ok for Multigroup_Frequency "
	  << "specialization.";
	PASSMSG(m.str());
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	gray_test();
	mg_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstOpacity_Builder_Helper, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstOpacity_Builder_Helper Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstOpacity_Builder_Helper on " << rtt_c4::node() 
	 << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstOpacity_Builder_Helper.cc
//---------------------------------------------------------------------------//
