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
#include "../Opacity.hh"
#include "../Mat_State.hh"
#include "../Global.hh"
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
using rtt_imc_test::IMC_Flat_Interface;
using rtt_imc::Flat_Mat_State_Builder;
using rtt_imc::Mat_State_Builder;
using rtt_imc::Mat_State;
using rtt_imc::Opacity;
using rtt_mc::OS_Builder;
using rtt_mc::OS_Mesh;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::OS_Mesh MT;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void flat_mat_state_test()
{
    // build a parser, mesh, and interface
    SP<Parser>             parser(new Parser("OS_Input"));
    SP<OS_Builder>         mb(new OS_Builder(parser));
    SP<OS_Mesh>            mesh = mb->build_Mesh();
    SP<IMC_Flat_Interface> interface(new IMC_Flat_Interface(mb));
    
    // pointer to a mat state builder
    SP<Mat_State_Builder<MT> > builder;

    // make a flat mat state builder
    builder = new Flat_Mat_State_Builder<MT>(interface);

    // check typeinfo
    if (typeid(*builder) != typeid(Flat_Mat_State_Builder<MT>))  ITFAILS;
    if (typeid(builder.bp()) != typeid(Mat_State_Builder<MT> *)) ITFAILS;

    // make a Mat_State
    SP<Mat_State<MT> > mat_state = builder->build_Mat_State(mesh);

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
	PASSMSG("Flat Mat_State passes all tests.");

    // make an Opacity
    SP<Opacity<MT> > opacity = builder->build_Opacity(mesh, mat_state);

    // check the opacity
    if (opacity->num_cells() != 6) ITFAILS;

    for (int cell = 1; cell <= 3; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double planck  = opacity->get_planck(cell);
	double fleck   = opacity->get_fleck(cell);
	double T       = mat_state->get_T(cell);
	double dedT    = mat_state->get_dedt(cell);
	double fplanck = opacity->fplanck(cell);
	double effabs  = opacity->get_sigeffabs(cell);
	double effsc   = opacity->get_sigeffscat(cell);
	
	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*1000.0*
	     mesh->volume(cell)*0.001/dedT*0.1);

	double fplanck_ref = fleck_ref * 0.1;
	double effabs_ref  = fleck_ref * 0.1;
	double effsc_ref   = (1.0 - fleck_ref) * 0.1;

	if (!soft_equiv(sig_abs, .1))          ITFAILS;
	if (!soft_equiv(sig_sc, .5))           ITFAILS;
	if (!soft_equiv(planck, .1))           ITFAILS;
	if (!soft_equiv(fleck, fleck_ref))     ITFAILS;
	if (!soft_equiv(fplanck, fplanck_ref)) ITFAILS;
	if (!soft_equiv(effabs, effabs_ref))   ITFAILS;
	if (!soft_equiv(effsc, effsc_ref))     ITFAILS;
    }

    for (int cell = 4; cell <= 6; cell++)
    {
	double sig_abs = opacity->get_sigma_abs(cell);
	double sig_sc  = opacity->get_sigma_thomson(cell);
	double planck  = opacity->get_planck(cell);
	double fleck   = opacity->get_fleck(cell);
	double T       = mat_state->get_T(cell);
	double dedT    = mat_state->get_dedt(cell);
	double fplanck = opacity->fplanck(cell);
	double effabs  = opacity->get_sigeffabs(cell);
	double effsc   = opacity->get_sigeffscat(cell);
	
	double fleck_ref   = 1.0 / 
	    (1.0 + 4.0*rtt_mc::global::a*rtt_mc::global::c*8000.0*
	     mesh->volume(cell)*0.001/dedT*0.02);

	double fplanck_ref = fleck_ref * 0.02;
	double effabs_ref  = fleck_ref * 0.02;
	double effsc_ref   = (1.0 - fleck_ref) * 0.02;

	if (!soft_equiv(sig_abs, .02))         ITFAILS;
	if (!soft_equiv(sig_sc, 0.0))          ITFAILS;
	if (!soft_equiv(planck, .02))          ITFAILS;
	if (!soft_equiv(fleck, fleck_ref))     ITFAILS;
	if (!soft_equiv(fplanck, fplanck_ref)) ITFAILS;
	if (!soft_equiv(effabs, effabs_ref))   ITFAILS;
	if (!soft_equiv(effsc, effsc_ref))     ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Flat Opacity passes all tests.");
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
	flat_mat_state_test();
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
