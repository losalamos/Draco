//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstCDI_Mat_State_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Mar  5 15:55:29 2003
 * \brief  test CDI_Mat_State_Builder class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Release.hh"
#include "../CDI_Mat_State_Builder.hh"
#include "../Opacity.hh"
#include "../Diffusion_Opacity.hh"
#include "../Mat_State.hh"
#include "../Frequency.hh"
#include "../Global.hh"
#include "../Particle.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
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

using namespace std;

using rtt_imc_test::Parser;
using rtt_imc_test::IMC_CDI_Interface;
using rtt_imc_test::IMC_CDI_Diffusion_Interface;
using rtt_imc::CDI_Mat_State_Builder;
using rtt_imc::Mat_State_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Diffusion_Opacity;
using rtt_imc::Opacity;
using rtt_mc::OS_Builder;
using rtt_mc::OS_Mesh;
using rtt_cdi::CDI;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::OS_Mesh               MT;
typedef rtt_imc::Gray_Frequency       G;
typedef rtt_imc::Multigroup_Frequency MG;
typedef rtt_imc::Particle<MT>         PT;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void gray_cdi_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                  parser(new Parser("OS_Input"));
    SP<OS_Builder>              mb(new OS_Builder(parser));
    SP<OS_Mesh>                 mesh = mb->build_Mesh();
    SP<IMC_CDI_Interface<PT> >  interface(new IMC_CDI_Interface<PT>());
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,G> > builder;

    // make a CDI mat state builder
    builder = new CDI_Mat_State_Builder<MT,G>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(CDI_Mat_State_Builder<MT,G>))   ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT,G> *)) ITFAILS;
    if (typeid(*builder) == typeid(CDI_Mat_State_Builder<MT,MG>))  ITFAILS;

    // build mat classes
    builder->build_mat_classes(mesh);

    // get the frequency
    SP<G> frequency = builder->get_Frequency();
    
    // brief check
    if (!frequency->is_gray()) ITFAILS;

    // get the Mat_State
    SP<Mat_State<MT> > mat_state = builder->get_Mat_State();

    // check the mat state
    {
	if (mat_state->num_cells() != 6) ITFAILS;

	if (!soft_equiv(mat_state->get_rho(1), 1.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(2), 2.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(3), 1.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(4), 3.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(5), 1.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(6), 2.0)) ITFAILS;

	for (int i = 1; i <= 6; i++)
	    if (!soft_equiv(mat_state->get_T(i), 3.0)) ITFAILS;

	if (!soft_equiv(mat_state->get_spec_heat(1), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(2), 2.7)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(3), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(4), 0.2)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(5), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(6), 0.2)) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("CDI Mat_State passes all tests in Gray test.");

    // get Opacity
    SP<Opacity<MT,G> > opacity = builder->get_Opacity();

    // check the opacities
    {
	if (opacity->num_cells() != 6) ITFAILS;
	
	double opac_1 = 100.0 / 27.0;
	double opac_2 = 1.5;
	double opac_3 = 1.0 + .1 * 3.0;
	double opac_s = 0.0;

	if (!soft_equiv(opacity->get_sigma_abs(1), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2), opac_3 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4), opac_2 * 3.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6), opac_2 * 2.0)) ITFAILS;

	for (int i = 1; i <= 6; i++)
	    if (!soft_equiv(opacity->get_sigma_thomson(i), opac_s)) ITFAILS;

	for (int i = 1; i <= 6; i++)
	{
	    double dedT = mat_state->get_dedt(i);
	    double opac = opacity->get_sigma_abs(i);
	    double fleck_ref   = 1.0 / 
		(1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*27.0*
		 mesh->volume(i)*0.001/dedT*opac);

	    if (!soft_equiv(opacity->get_fleck(i), fleck_ref)) ITFAILS;

	    double effs = (1.0 - fleck_ref) * opac;
	    double effa = fleck_ref * opac;

	    if (!soft_equiv(opacity->get_sigeffscat(i), effs)) ITFAILS;
	    if (!soft_equiv(opacity->get_sigeffabs(i), effa))  ITFAILS;
	}

	for (int i = 1; i <= 6; i++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(i),1.0))
		ITFAILS;
    }

    // diffusion opacity should be null
    SP<Diffusion_Opacity<MT> > diff = builder->get_Diffusion_Opacity();
    if (diff) ITFAILS;

    if (rtt_imc_test::passed)
	PASSMSG("CDI Opacity passes all tests in Gray test.");
}
//---------------------------------------------------------------------------//

void gray_cdi_diffusion_mat_state_test()
{
    // do test with analytic data defined in CDI

    // build a parser, mesh, and interface
    SP<Parser>                  parser(new Parser("OS_Input"));
    SP<OS_Builder>              mb(new OS_Builder(parser));
    SP<OS_Mesh>                 mesh = mb->build_Mesh();
    SP<IMC_CDI_Interface<PT> >  interface(new IMC_CDI_Interface<PT>(1));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,G> > builder;

    // make a CDI mat state builder
    builder = new CDI_Mat_State_Builder<MT,G>(interface);

    // build mat classes
    builder->build_mat_classes(mesh);

    // get opacity and mat state
    SP<Mat_State<MT> >  mat_state = builder->get_Mat_State();
    SP<Opacity<MT, G> > opacity   = builder->get_Opacity();

    // get the diffusion opacities
    SP<Diffusion_Opacity<MT> > diff = builder->get_Diffusion_Opacity();
    if (!diff)                  ITFAILS;
    if (diff->num_cells() != 6) ITFAILS;

    // check the opacities, the rosseland opacities should be the absorption
    // + scattering because these opacities are analytic
    {
	double opac_s = 1.01;
	double opac_1 = 100.0 / 27.0 + opac_s;
	double opac_2 = 1.5 + opac_s;
	double opac_3 = 1.0 + .1 * 3.0 + opac_s;

	if (!soft_equiv(diff->get_Rosseland_opacity(1), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(2), opac_3 * 2.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(3), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(4), opac_2 * 3.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(5), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(6), opac_2 * 2.0)) ITFAILS;

	for (int i = 1; i <= 6; i++)
	{
	    double dedT = mat_state->get_dedt(i);
	    double opac = opacity->get_sigma_abs(i);
	    double fleck_ref   = 1.0 / 
		(1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*27.0*
		 mesh->volume(i)*0.001/dedT*opac);

	    if (!soft_equiv(diff->get_fleck(i), fleck_ref)) ITFAILS;
	}
    }

    if (rtt_imc_test::passed)
    {
	ostringstream m;
	m << "CDI Diffusion_Opacity passes tests for interface with"
	  << " analytic CDI opacity.";
	PASSMSG(m.str());
    }

    // do test with rosseland opacities defined in the CDI

    // define new interface
    SP<IMC_CDI_Diffusion_Interface<PT> > diff_interface(
	new IMC_CDI_Diffusion_Interface<PT>("gray"));

    // make a new builder
    builder = new CDI_Mat_State_Builder<MT,G>(diff_interface);

    // build mat classes
    builder->build_mat_classes(mesh);

    // get opacity and mat state
    mat_state = builder->get_Mat_State();
    opacity   = builder->get_Opacity();

    // get the diffusion opacities
    diff = builder->get_Diffusion_Opacity();
    if (!diff)                  ITFAILS;
    if (diff->num_cells() != 6) ITFAILS;

    // check the opacities, the rosseland opacity is part of the CDI and has
    // a value of 1.1 + 1/T^3
    {
	double opac_1 = 1.1 + 1.0 / 27.0;

	if (!soft_equiv(diff->get_Rosseland_opacity(1), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(2), opac_1 * 2.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(3), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(4), opac_1 * 3.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(5), opac_1 * 1.0)) ITFAILS;
	if (!soft_equiv(diff->get_Rosseland_opacity(6), opac_1 * 2.0)) ITFAILS;

	for (int i = 1; i <= 6; i++)
	{
	    double dedT = mat_state->get_dedt(i);
	    double opac = opacity->get_sigma_abs(i);
	    double fleck_ref   = 1.0 / 
		(1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*27.0*
		 mesh->volume(i)*0.001/dedT*opac);

	    if (!soft_equiv(diff->get_fleck(i), fleck_ref)) ITFAILS;
	}
    }

    if (rtt_imc_test::passed)
    {
	ostringstream m;
	m << "CDI Diffusion_Opacity passes tests for interface with"
	  << " Rosseland CDI opacities defined.";
	PASSMSG(m.str());
    }

    // now check that if the only gray opacities we have in CDI are Planck
    // the problem should throw an assertion
    SP<IMC_CDI_Diffusion_Interface<PT> > bad_interface(
	new IMC_CDI_Diffusion_Interface<PT>("gray", true));
    bool caught = false;
    try
    {
	// make a new builder
	builder = new CDI_Mat_State_Builder<MT,G>(bad_interface);
	
	// build mat classes
	builder->build_mat_classes(mesh);	
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream m;
	m << "Good, caught the following assertion, " << ass.what();
	PASSMSG(m.str());
    }
    if (!caught)
    {
	FAILMSG("Failed to catch assertion in CDI_Mat_State_Builder.");
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("CDI Diffusion_Opacity passes all tests in Gray test.");
}

//---------------------------------------------------------------------------//

void mg_cdi_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                   parser(new Parser("OS_Input"));
    SP<OS_Builder>               mb(new OS_Builder(parser));
    SP<OS_Mesh>                  mesh = mb->build_Mesh();
    SP<IMC_CDI_Interface<PT> >   interface(new IMC_CDI_Interface<PT>());
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,MG> > builder;

    // make a CDI mat state builder
    builder = new CDI_Mat_State_Builder<MT,MG>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(CDI_Mat_State_Builder<MT,MG>))   ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT,MG> *)) ITFAILS;
    if (typeid(*builder) == typeid(CDI_Mat_State_Builder<MT,G>))    ITFAILS;

    // build mat classes
    builder->build_mat_classes(mesh);

    // build the frequency
    SP<MG> frequency = builder->get_Frequency();
    
    // brief check
    if (!frequency->is_multigroup()) ITFAILS;

    // make a Mat_State
    SP<Mat_State<MT> > mat_state = builder->get_Mat_State();

    // check the mat state
    {
	if (mat_state->num_cells() != 6) ITFAILS;

	if (!soft_equiv(mat_state->get_rho(1), 1.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(2), 2.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(3), 1.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(4), 3.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(5), 1.0)) ITFAILS;
	if (!soft_equiv(mat_state->get_rho(6), 2.0)) ITFAILS;

	for (int i = 1; i <= 6; i++)
	    if (!soft_equiv(mat_state->get_T(i), 3.0)) ITFAILS;

	if (!soft_equiv(mat_state->get_spec_heat(1), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(2), 2.7)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(3), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(4), 0.2)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(5), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(6), 0.2)) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("CDI Mat_State passes all tests in MG test.");

    // make an Opacity
    SP<Opacity<MT,MG> > opacity = builder->get_Opacity();

    // check the groups
    vector<double> ref_groups(4);
    ref_groups[0] = 0.01;
    ref_groups[1] = 0.1;
    ref_groups[2] = 1.0;
    ref_groups[3] = 10.0;
    {
	vector<double> groups = frequency->get_group_boundaries();
	if (!soft_equiv(groups.begin(), groups.end(), 
			ref_groups.begin(), ref_groups.end())) ITFAILS;
    }

    // check the opacity
    {
	if (!soft_equiv(opacity->get_sigma_abs(1, 1), 5.0 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 2), 1.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 3), 0.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 1), 3.0 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 2), 1.3 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 3), 0.3 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 1), 5.0 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 2), 1.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 3), 0.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 1), 4.0 * 3.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 2), 1.4 * 3.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 3), 0.4 * 3.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 1), 5.0 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 2), 1.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 3), 0.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 1), 4.0 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 2), 1.4 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 3), 0.4 * 2.0)) ITFAILS;

	for (int c = 1; c <= 6; c++)
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
		if (!soft_equiv(opacity->get_sigma_thomson(c,g), 0.0)) 
		    ITFAILS;
	
	// check integrated norm Planck
	double int_Planck = CDI::integratePlanckSpectrum(0.01, 10.0, 3.0);
	
	for (int c = 1; c <= 6; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	// check effective cross sections
	for (int c = 1; c < 6; c++)
	{
	    pair<double,double> b;
	    double              sum = 0.0;
	    vector<double>      ref_emission(3);
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		b                 = frequency->get_group_boundaries(g);
		sum              += opacity->get_sigma_abs(c,g) *
		    CDI::integratePlanckSpectrum(b.first, b.second, 3.0);
		ref_emission[g-1] = sum;
	    }
	    
	    vector<double> emission = opacity->get_emission_group_cdf(c);

	    if (!soft_equiv(emission.begin(), emission.end(),
			    ref_emission.begin(), ref_emission.end())) ITFAILS;

	    double planck = sum / opacity->get_integrated_norm_Planck(c);

	    double beta   = 4.0 * rtt_mc::global::a * 27.0 * mesh->volume(c) /
		mat_state->get_dedt(c);

	    double fleck  = 1.0 / 
		(1.0 + .001 * rtt_mc::global::c * beta * planck);  

	    if (!soft_equiv(opacity->get_fleck(c), fleck)) ITFAILS;

	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		double effabs = opacity->get_sigma_abs(c,g) * fleck;
		double effsct = opacity->get_sigma_abs(c,g) * (1.0-fleck);

		if (!soft_equiv(opacity->get_sigeffscat(c,g),effsct)) ITFAILS;
		if (!soft_equiv(opacity->get_sigeffabs(c,g),effabs))  ITFAILS;
	    }
	}
    }

    if (rtt_imc_test::passed)
	PASSMSG("CDI MG Frequency and MG Opacity passes all tests.");
}

//---------------------------------------------------------------------------//

void mg_cdi_diffusion_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                   parser(new Parser("OS_Input"));
    SP<OS_Builder>               mb(new OS_Builder(parser));
    SP<OS_Mesh>                  mesh = mb->build_Mesh();

    // build dummy interface with hybrid diffusion on and common mg opacities
    // off
    SP<IMC_CDI_Interface<PT> >   interface(new IMC_CDI_Interface<PT>(1));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,MG> > builder;

    // make a CDI mat state builder
    builder = new CDI_Mat_State_Builder<MT,MG>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(CDI_Mat_State_Builder<MT,MG>))   ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT,MG> *)) ITFAILS;
    if (typeid(*builder) == typeid(CDI_Mat_State_Builder<MT,G>))    ITFAILS;

    // build mat classes
    builder->build_mat_classes(mesh);

    // build the mat objects
    SP<MG>                     frequency = builder->get_Frequency();
    SP<Mat_State<MT> >         mat_state = builder->get_Mat_State();
    SP<Opacity<MT,MG> >        opacity   = builder->get_Opacity();
    SP<Diffusion_Opacity<MT> > diff      = builder->get_Diffusion_Opacity();
    if (!diff)                  ITFAILS;
    if (diff->num_cells() != 6) ITFAILS;

    // check the groups
    vector<double> ref_groups(4);
    ref_groups[0] = 0.01;
    ref_groups[1] = 0.1;
    ref_groups[2] = 1.0;
    ref_groups[3] = 10.0;
    {
	vector<double> groups = frequency->get_group_boundaries();
	if (!soft_equiv(groups.begin(), groups.end(), 
			ref_groups.begin(), ref_groups.end())) ITFAILS;
    }

    // check the opacity
    double int_Planck;
    {
	if (!soft_equiv(opacity->get_sigma_abs(1, 1), 5.0 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 2), 1.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 3), 0.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 1), 3.0 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 2), 1.3 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 3), 0.3 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 1), 5.0 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 2), 1.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 3), 0.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 1), 4.0 * 3.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 2), 1.4 * 3.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 3), 0.4 * 3.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 1), 5.0 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 2), 1.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 3), 0.5 * 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 1), 4.0 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 2), 1.4 * 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 3), 0.4 * 2.0)) ITFAILS;

	for (int c = 1; c <= 6; c++)
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
		if (!soft_equiv(opacity->get_sigma_thomson(c,g), 
				1.01*mat_state->get_rho(c))) ITFAILS;
	
	// check integrated norm Planck
	int_Planck = CDI::integratePlanckSpectrum(0.01, 10.0, 3.0);
	for (int c = 1; c <= 6; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	// check effective cross sections in opacity and diffusion opacity
	for (int c = 1; c < 6; c++)
	{
	    pair<double,double> b;
	    double              sum = 0.0;
	    vector<double>      ref_emission(3);
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		b                 = frequency->get_group_boundaries(g);
		sum              += opacity->get_sigma_abs(c,g) *
		    CDI::integratePlanckSpectrum(b.first, b.second, 3.0);
		ref_emission[g-1] = sum;
	    }
	    
	    vector<double> emission = opacity->get_emission_group_cdf(c);

	    if (!soft_equiv(emission.begin(), emission.end(),
			    ref_emission.begin(), ref_emission.end())) ITFAILS;

	    double planck = sum / opacity->get_integrated_norm_Planck(c);

	    double beta   = 4.0 * rtt_mc::global::a * 27.0 * mesh->volume(c) /
		mat_state->get_dedt(c);

	    double fleck  = 1.0 / 
		(1.0 + .001 * rtt_mc::global::c * beta * planck);  

	    if (!soft_equiv(opacity->get_fleck(c), fleck)) ITFAILS;

	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		double effabs = opacity->get_sigma_abs(c,g) * fleck;
		double effsct = opacity->get_sigma_abs(c,g) * (1.0-fleck);

		if (!soft_equiv(opacity->get_sigeffscat(c,g),effsct)) ITFAILS;
		if (!soft_equiv(opacity->get_sigeffabs(c,g),effabs))  ITFAILS;
	    }

	    // check Rosseland opacities in diffusion opacity
	    double ros_sum         = 0.0;
	    double inv_sig_ros_sum = 0.0;
	    double r_g             = 0.0;
	    double ros_ref         = 0.0;
	    double T               = 3.0;
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		b   = frequency->get_group_boundaries(g);
		r_g = CDI::integrateRosselandSpectrum(b.first, b.second, T);

		// sums (scattering is 1.01)
		ros_sum         += r_g;
		inv_sig_ros_sum += r_g / (opacity->get_sigma_abs(c, g) + 
					  1.01 * mat_state->get_rho(c));
	    }

	    // check ros sum
	    if (!soft_equiv(
		    CDI::integrateRosselandSpectrum(0.01, 10.0, T),
		    ros_sum)) ITFAILS;

	    // calculate reference rosseland and compare
	    ros_ref = ros_sum / inv_sig_ros_sum;

	    // check it
	    if (!soft_equiv(diff->get_Rosseland_opacity(c), ros_ref)) ITFAILS;
	}
    }

    if (rtt_imc_test::passed)
	PASSMSG("CDI Diffusion_Opacity passes all tests in MG test.");
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
	
	// gray tests
	gray_cdi_mat_state_test();
	gray_cdi_diffusion_mat_state_test();

	// mg tests
	mg_cdi_mat_state_test();
	mg_cdi_diffusion_mat_state_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstCDI_Mat_State_Builder, " << ass.what()
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
	    cout << "**** tstCDI_Mat_State_Builder Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstCDI_Mat_State_Builder on " 
	 << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstCDI_Mat_State_Builder.cc
//---------------------------------------------------------------------------//
