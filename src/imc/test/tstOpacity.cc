//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstOpacity.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 19 11:23:44 2003
 * \brief  Test imc-based opacity classes.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Release.hh"
#include "../Diffusion_Opacity.hh"
#include "../Opacity.hh"
#include "../Fleck_Factors.hh"
#include "../Frequency.hh"
#include "mc/Constants.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_imc_test::Parser;
using rtt_imc::Diffusion_Opacity;
using rtt_imc::Opacity;
using rtt_imc::Fleck_Factors;
using rtt_imc::Gray_Frequency;
using rtt_imc::Multigroup_Frequency;
using rtt_mc::OS_Builder;
using rtt_mc::OS_Mesh;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::OS_Mesh           MT;
typedef MT::CCSF<double>          ccsf_double;
typedef MT::CCSF<vector<double> > ccsf_sf_double;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

SP<MT> build_mesh()
{
    // build a parser, mesh, and interface
    SP<Parser>      parser(new Parser("OS_Input"));
    SP<OS_Builder>  mb(new OS_Builder(parser));
    SP<OS_Mesh>     mesh = mb->build_Mesh();

    if (mesh->num_cells() != 6)             ITFAILS;
    if (mesh->get_spatial_dimension() != 2) ITFAILS;

    if (rtt_imc_test::passed)
	PASSMSG("Built simple 2D mesh for testing opacities.");

    return mesh;
}

//---------------------------------------------------------------------------//

void test_diffusion_opacity()
{
    // get the mesh
    SP<MT> mesh = build_mesh();

    // make a rosseland field
    ccsf_double r(mesh);

    // make a rosseland scattering field
    ccsf_double s(mesh);
    
    double factor = 10.5;
    for (ccsf_double::iterator itr = r.begin(); itr != r.end(); itr++)
    {
	*itr    = 1.1 * factor;
	factor *= 2.0;
    }
    factor = 10.5;

    for (ccsf_double::iterator itr = s.begin(); itr != s.end(); itr++)
	*itr = 1.01;

    // make fleck factors
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    for (int i = 1; i <= fleck->fleck.size(); i++)
	fleck->fleck(i) = 0.023 + 0.1 * static_cast<double>(i);

    // make a diffusion opacity
    Diffusion_Opacity<MT> opac(fleck, r, s);

    // now check
    if (opac.num_cells() != 6) ITFAILS;

    for (int i = 1; i <= 6; i++)
    {
	double f = 0.023 + 0.1 * static_cast<double>(i);
	double r = 1.1 * factor;
	double s = 1.01;
	double D = rtt_mc::global::c / (3.0 * ((1.0-f) * r + s));
	double e = 1.0 / ((1.0 - f) * r + s);

	if (!soft_equiv(opac.get_fleck(i), f))             ITFAILS;
	if (!soft_equiv(opac.get_Rosseland_opacity(i), r)) ITFAILS;
	if (!soft_equiv(opac.get_gray_scattering(i), s))   ITFAILS;
	if (!soft_equiv(opac.get_random_walk_D(i), D))     ITFAILS;
	if (!soft_equiv(opac.get_Rosseland_effmfp(i), e))  ITFAILS;

	factor  *= 2.0;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Diffusion_Opacity ok.");
}

//---------------------------------------------------------------------------//

void test_gray_opacity()
{
    // get the mesh
    SP<MT> mesh = build_mesh();

    // make frequency
    SP<Gray_Frequency> gray(new Gray_Frequency);

    // make opacity fields
    ccsf_double abs(mesh);
    ccsf_double sct(mesh);

    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	abs(i) = 0.125;
	sct(i) = 0.1;
    }

    // make fleck factors
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    for (int i = 1; i <= fleck->fleck.size(); i++)
	fleck->fleck(i) = 0.023 + 0.1 * static_cast<double>(i);

    // build opacity
    Opacity<MT, Gray_Frequency> opac(gray, abs, sct, fleck);

    // testit
    if (opac.num_cells() != 6) ITFAILS;

    for (int i = 1; i <= 6; i++)
    {
	double f  = 0.023 + 0.1 * static_cast<double>(i); 
	double ea = f * 0.125;
	double es = (1.0-f) * 0.125;

	if (!soft_equiv(opac.get_fleck(i), f))           ITFAILS;
	if (!soft_equiv(opac.get_sigma_abs(i), 0.125))   ITFAILS;
	if (!soft_equiv(opac.get_sigma_thomson(i), 0.1)) ITFAILS;
	if (!soft_equiv(opac.get_sigeffscat(i), es))     ITFAILS;
	if (!soft_equiv(opac.get_sigeffabs(i), ea))      ITFAILS;

	if (!opac.get_Frequency()->is_gray())            ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Gray Opacity specialization ok.");
}

//---------------------------------------------------------------------------//

void test_mg_opacity()
{
    // get the mesh
    SP<MT> mesh = build_mesh();

    // group bounds
    vector<double> g(4, 0.0);
    {
	g[0] = 0.001;
	g[1] = 0.01;
	g[2] = 1.0;
	g[3] = 10.0;
    }
    SP<Multigroup_Frequency> mg(new Multigroup_Frequency(g));

    if (mg->get_num_groups() != 3) ITFAILS;

    // make fleck factors
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    for (int i = 1; i <= fleck->fleck.size(); i++)
	fleck->fleck(i) = 0.023 + 0.1 * static_cast<double>(i);

    // opacities
    ccsf_sf_double abs(mesh);
    ccsf_sf_double sct(mesh);
    
    // integrators
    ccsf_double    planck(mesh);
    ccsf_sf_double cdf(mesh);

    double cf = 0.1;
    for (int i = 1; i <= 6; i++)
    {
	abs(i).resize(3);
	sct(i).resize(3);
	cdf(i).resize(3);

	for (int j = 0; j < 3; j++)
	{
	    abs(i)[j] = 0.125 + cf;
	    sct(i)[j] = 0.1 + cf;
	    cdf(i)[j] = 0.2 + cf;
	}

	planck(i) = 0.1 + cf;

	cf += 0.1;
    }

    Opacity<MT, Multigroup_Frequency> opac(mg, abs, sct, fleck, planck, cdf);

    // testit
    if (opac.num_cells() != 6) ITFAILS;

    cf = 0.1;
    for (int i = 1; i <= 6; i++)
    {
	double f = 0.023 + 0.1 * static_cast<double>(i); 
	double p = 0.1 + cf; 

	if (!soft_equiv(opac.get_fleck(i), f))                  ITFAILS;
	if (!soft_equiv(opac.get_integrated_norm_Planck(i), p)) ITFAILS; 

	for (int g = 1; g <= 3; g++)
	{
	    double ea  = f * (0.125 + cf);
	    double es  = (1.0-f) * (0.125 + cf);
	    
	    if (!soft_equiv(opac.get_sigma_abs(i,g), 0.125+cf))   ITFAILS;
	    if (!soft_equiv(opac.get_sigma_thomson(i,g), 0.1+cf)) ITFAILS;
	    if (!soft_equiv(opac.get_sigeffscat(i,g), es))        ITFAILS;
	    if (!soft_equiv(opac.get_sigeffabs(i,g), ea))         ITFAILS;
	}

        double cdf = 0.2 + cf;
	if (!soft_equiv(opac.get_integrated_sigma_times_Planck(i), cdf)) 
	    ITFAILS;

	if (!opac.get_Frequency()->is_multigroup())     ITFAILS;

	if (opac.get_emission_group_cdf(i).size() != 3) ITFAILS;

	cf += 0.1;
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("Multigroup Opacity specialization ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // this is an inherently serial test
    if (rtt_c4::node())
    {
	rtt_c4::finalize();
	return 0;
    }

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
	test_diffusion_opacity();
	test_gray_opacity();
	test_mg_opacity();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstOpacity, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstOpacity Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstOpacity on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstOpacity.cc
//---------------------------------------------------------------------------//
