//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstMat_State_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 16:51:58 2001
 * \brief  Mat_State_Builder test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Release.hh"
#include "../Flat_Mat_State_Builder.hh"
#include "../CDI_Mat_State_Builder.hh"
#include "../Opacity.hh"
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

using namespace std;

using rtt_imc_test::Parser;
using rtt_imc_test::IMC_Flat_Interface;
using rtt_imc_test::IMC_CDI_Interface;
using rtt_imc::Flat_Mat_State_Builder;
using rtt_imc::CDI_Mat_State_Builder;
using rtt_imc::Mat_State_Builder;
using rtt_imc::Mat_State;
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

void gray_flat_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                  parser(new Parser("OS_Input"));
    SP<OS_Builder>              mb(new OS_Builder(parser));
    SP<OS_Mesh>                 mesh = mb->build_Mesh();
    SP<IMC_Flat_Interface<PT> > interface(new IMC_Flat_Interface<PT>(mb));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,G> > builder;

    // make a flat mat state builder
    builder = new Flat_Mat_State_Builder<MT,G>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(Flat_Mat_State_Builder<MT,G>))  ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT,G> *)) ITFAILS;
    if (typeid(*builder) == typeid(Flat_Mat_State_Builder<MT,MG>)) ITFAILS;

    // build objects
    builder->build_mat_classes(mesh);

    // build the frequency
    SP<G> frequency = builder->get_Frequency();
    
    // brief check
    if (!frequency->is_gray()) ITFAILS;

    // make a Mat_State
    SP<Mat_State<MT> > mat_state = builder->get_Mat_State();

    // check the mat_state
    if (mat_state->num_cells() != 6) ITFAILS;

    for (int cell = 1; cell <= 3; cell++)
    {
	double d        = mat_state->get_rho(cell);
	double T        = mat_state->get_T(cell);
	double Cv       = mat_state->get_spec_heat(cell);
	double dedT     = mat_state->get_dedt(cell);
	double dedT_ref = 0.1 * 1.0 * mesh->volume(cell);

	if (!soft_equiv(d, 1.0))         ITFAILS;
	if (!soft_equiv(T, 10.0))        ITFAILS;
	if (!soft_equiv(Cv, 0.1))        ITFAILS;
	if (!soft_equiv(dedT, dedT_ref)) ITFAILS;
    }
    for (int cell = 4; cell <= 6; cell++)
    {
	double d        = mat_state->get_rho(cell);
	double T        = mat_state->get_T(cell);
	double Cv       = mat_state->get_spec_heat(cell);
	double dedT     = mat_state->get_dedt(cell);
	double dedT_ref = 0.2 * 2.0 * mesh->volume(cell);

	if (!soft_equiv(d, 2.0))         ITFAILS;
	if (!soft_equiv(T, 20.0))        ITFAILS;
	if (!soft_equiv(Cv, 0.2))        ITFAILS;
	if (!soft_equiv(dedT, dedT_ref)) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat Mat_State passes all tests in Gray test.");

    // make an Opacity
    SP<Opacity<MT,G> > opacity = builder->get_Opacity();

    // check the opacity
    if (opacity->num_cells() != 6) ITFAILS;

    for (int cell = 1; cell <= 3; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double fleck   = opacity->get_fleck(cell);
	double T       = mat_state->get_T(cell);
	double dedT    = mat_state->get_dedt(cell);
	double effabs  = opacity->get_sigeffabs(cell);
	double effsc   = opacity->get_sigeffscat(cell);
	
	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*1000.0*
	     mesh->volume(cell)*0.001/dedT*0.1);

	double effabs_ref  = fleck_ref * 0.1;
	double effsc_ref   = (1.0 - fleck_ref) * 0.1;

	if (!soft_equiv(sig_abs, .1))          ITFAILS;
	if (!soft_equiv(sig_sc, .5))           ITFAILS;
	if (!soft_equiv(fleck, fleck_ref))     ITFAILS;
	if (!soft_equiv(effabs, effabs_ref))   ITFAILS;
	if (!soft_equiv(effsc, effsc_ref))     ITFAILS;
    }

    for (int cell = 4; cell <= 6; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double fleck   = opacity->get_fleck(cell);
	double T       = mat_state->get_T(cell);
	double dedT    = mat_state->get_dedt(cell);
	double effabs  = opacity->get_sigeffabs(cell);
	double effsc   = opacity->get_sigeffscat(cell);
	
	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*8000.0*
	     mesh->volume(cell)*0.001/dedT*0.02);

	double effabs_ref  = fleck_ref * 0.02;
	double effsc_ref   = (1.0 - fleck_ref) * 0.02;

	if (!soft_equiv(sig_abs, .02))         ITFAILS;
	if (!soft_equiv(sig_sc, 0.0))          ITFAILS;
	if (!soft_equiv(fleck, fleck_ref))     ITFAILS;
	if (!soft_equiv(effabs, effabs_ref))   ITFAILS;
	if (!soft_equiv(effsc, effsc_ref))     ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat Opacity passes all tests in Gray test.");
}

//---------------------------------------------------------------------------//

void mg_flat_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>                    parser(new Parser("OS_Input"));
    SP<OS_Builder>                mb(new OS_Builder(parser));
    SP<OS_Mesh>                   mesh = mb->build_Mesh();
    SP<IMC_Flat_Interface<PT> >   interface(new IMC_Flat_Interface<PT>(mb));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT,MG> > builder;

    // make a CDI mat state builder
    builder = new Flat_Mat_State_Builder<MT,MG>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(Flat_Mat_State_Builder<MT,MG>))  ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT,MG> *)) ITFAILS;
    if (typeid(*builder) == typeid(Flat_Mat_State_Builder<MT,G>))   ITFAILS;

    // build objects
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

	if (!soft_equiv(mat_state->get_rho(1), 1.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(2), 1.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(3), 1.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(4), 2.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(5), 2.0))       ITFAILS;
	if (!soft_equiv(mat_state->get_rho(6), 2.0))       ITFAILS;

	if (!soft_equiv(mat_state->get_T(1), 10.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(2), 10.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(3), 10.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(4), 20.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(5), 20.0))        ITFAILS;
	if (!soft_equiv(mat_state->get_T(6), 20.0))        ITFAILS;

	if (!soft_equiv(mat_state->get_spec_heat(1), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(2), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(3), 0.1)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(4), 0.2)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(5), 0.2)) ITFAILS;
	if (!soft_equiv(mat_state->get_spec_heat(6), 0.2)) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat Mat_State passes all tests in MG test.");

    // make an Opacity
    SP<Opacity<MT,MG> > opacity = builder->get_Opacity();

    // check the groups
    vector<double> ref_groups(4);
    ref_groups[0] = 0.01;
    ref_groups[1] = 0.1;
    ref_groups[2] = 15.0;
    ref_groups[3] = 100.0;
    {
	vector<double> groups = frequency->get_group_boundaries();
	if (!soft_equiv(groups.begin(), groups.end(), 
			ref_groups.begin(), ref_groups.end())) ITFAILS;
    }

    // check the opacity
    {
	if (!soft_equiv(opacity->get_sigma_abs(1, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(1, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(2, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 1), 1.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 2), 0.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(3, 3), 0.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(4, 3), 1.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(5, 3), 1.1)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 1), 2.0)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 2), 1.5)) ITFAILS;
	if (!soft_equiv(opacity->get_sigma_abs(6, 3), 1.1)) ITFAILS;

	for (int c = 1; c <= 6; c++)
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
		if (!soft_equiv(opacity->get_sigma_thomson(c,g), 0.0)) 
		    ITFAILS;
	
	// check integrated norm Planck
	double int_Planck;

	int_Planck = CDI::integratePlanckSpectrum(0.01, 100.0, 10.0);
	for (int c = 1; c <= 3; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	int_Planck = CDI::integratePlanckSpectrum(0.01, 100.0, 20.0);
	for (int c = 4; c <= 6; c++)
	    if (!soft_equiv(opacity->get_integrated_norm_Planck(c),
			    int_Planck)) ITFAILS;

	// check effective cross sections
	for (int c = 1; c < 6; c++)
	{
	    pair<double,double> b;
	    double              sum = 0.0;
	    vector<double>      ref_emission(3);
	    double              T = mat_state->get_T(c);
	    for (int g = 1; g <= frequency->get_num_groups(); g++)
	    {
		b                 = frequency->get_group_boundaries(g);
		sum              += opacity->get_sigma_abs(c,g) *
		    CDI::integratePlanckSpectrum(b.first, b.second, T);
		ref_emission[g-1] = sum;
	    }
	    
	    vector<double> emission = opacity->get_emission_group_cdf(c);

	    if (!soft_equiv(emission.begin(), emission.end(),
			    ref_emission.begin(), ref_emission.end())) ITFAILS;

	    double planck = sum / opacity->get_integrated_norm_Planck(c);

	    double beta   = 4.0 * rtt_mc::global::a * T*T*T * mesh->volume(c) /
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
	PASSMSG("Flat MG Frequency and MG Opacity passes all tests.");
}

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

    if (rtt_imc_test::passed)
	PASSMSG("CDI Opacity passes all tests in Gray test.");
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
	gray_flat_mat_state_test();
	mg_flat_mat_state_test();

	gray_cdi_mat_state_test();
	mg_cdi_mat_state_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstMat_State_Builder, " << ass.what()
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
	    cout << "**** tstMat_State_Builder Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstMat_State_Builder on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstMat_State_Builder.cc
//---------------------------------------------------------------------------//
