//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstOpacity.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:05:50 2001
 * \brief  Opacity and Mat_State tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Opacity_Builder.hh"
#include "../Opacity.hh"
#include "../Mat_State.hh"
#include "../Release.hh"
#include "../Global.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::global::soft_equiv;
using rtt_imc::Opacity;
using rtt_imc::Opacity_Builder;
using rtt_imc::Mat_State;
using rtt_imc_test::IMC_Interface;
using rtt_imc_test::Parser;
using rtt_dsxx::SP;

// some typedefs
typedef Mat_State<OS_Mesh> MSOS;
typedef Opacity<OS_Mesh> OOS;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// test Opacity Builder and Mat_State
void Mat_Test()
{
    // build a mesh
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // get an interface (dummy)
    SP<IMC_Interface> interface(new IMC_Interface(mb));

    // build the Mat_State
    Opacity_Builder<OS_Mesh> ob(interface);
    SP<MSOS> mat = ob.build_Mat(mesh);
    MSOS &m = *mat;

    // test Mat_State
    for (int i = 1; i <= 3; i++)
    {
	// density
	if (m.get_rho(i) != 1.0)   ITFAILS;
	if (m.get_rho(i+3) != 2.0) ITFAILS;

	// Temperature
	if (m.get_T(i) != 10.0)   ITFAILS;
	if (m.get_T(i+3) != 20.0) ITFAILS;

	// dedt = Cv * rho * V
	if (m.get_dedt(i) != .1 * 1.0 * 2.0)   ITFAILS;
	if (m.get_dedt(i+3) != .2 * 2.0 * 2.0) ITFAILS;
	
	// spec heat (Cv)
	if (m.get_spec_heat(i) != .1)   ITFAILS;
	if (m.get_spec_heat(i+3) != .2) ITFAILS;
    }
    
    if (m.num_cells() != 6) ITFAILS;
}

//---------------------------------------------------------------------------//
// test Opacity Builder, Mat_State, and Opacity
void Opacity_Test()
{
    // build a mesh
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // get an interface (dummy)
    SP<IMC_Interface> interface(new IMC_Interface(mb));

    // build the Mat_State and Opacity
    Opacity_Builder<OS_Mesh> ob(interface);
    SP<MSOS> mat    = ob.build_Mat(mesh);
    SP<OOS> opacity = ob.build_Opacity(mesh, mat);
    MSOS &m = *mat; 
    OOS &o  = *opacity;

    // test opacity
    for (int i = 1; i <= 3; i++)
    {
	// test abs opacity
	if (o.get_sigma_abs(i) != .1 * m.get_rho(i))      ITFAILS;
	if (o.get_sigma_abs(i+3) != .01 * m.get_rho(i+3)) ITFAILS;

	// test sigma thomson
	if (o.get_sigma_thomson(i) != .5 * m.get_rho(i))    ITFAILS;
	if (o.get_sigma_thomson(i+3) != 0 * m.get_rho(i+3)) ITFAILS;

	// test planck
	if (o.get_planck(i) != .1 * m.get_rho(i))      ITFAILS;
	if (o.get_planck(i+3) != .01 * m.get_rho(i+3)) ITFAILS; 
	
	// test fleck
	
	// calculate beta at i and beta and i+3
	double bi = 4.0 * rtt_mc::global::a * m.get_T(i) * m.get_T(i) *
	    m.get_T(i) * mesh->volume(i) / m.get_dedt(i);
	double bi3 = 4.0 * rtt_mc::global::a * m.get_T(i+3) * m.get_T(i+3) *
	    m.get_T(i+3) * mesh->volume(i+3) / m.get_dedt(i+3);
	double di = 1.0 * bi * rtt_mc::global::c * .001 
	    * o.get_sigma_abs(i);
	double di3 = 1.0 * bi3 * rtt_mc::global::c * .001
	    * o.get_sigma_abs(i+3);
	double flecki  = 1.0 / (1.0 + di);
	double flecki3 = 1.0 / (1.0 + di3);
	if (!soft_equiv(o.get_fleck(i), flecki))    ITFAILS;
	if (!soft_equiv(o.get_fleck(i+3), flecki3)) ITFAILS;

	// test fleck * planck -Wold-style-cast
	if (!soft_equiv(o.fplanck(i), flecki * o.get_planck(i), 1.e-16))
	    ITFAILS;
	if (!soft_equiv(o.fplanck(i+3), flecki3 * o.get_planck(i+3), 1.e-16)) 
	    ITFAILS;

	// test Fleck effective cross sections
	double escati  = (1 - flecki) * .1 * m.get_rho(i);
	double eabsi   = flecki * .1 * m.get_rho(i);
	double escati3 = (1 - flecki3) * .01 * m.get_rho(i+3);
	double eabsi3  = flecki3 * .01 * m.get_rho(i+3);
	if (!soft_equiv(o.get_sigeffscat(i), escati, 1.e-16))    ITFAILS;
	if (!soft_equiv(o.get_sigeffabs(i), eabsi, 1.e-16))      ITFAILS;
	if (!soft_equiv(o.get_sigeffscat(i+3), escati3, 1.e-16)) ITFAILS;
	if (!soft_equiv(o.get_sigeffabs(i+3), eabsi3, 1.e-16))   ITFAILS; 
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is essentially a scalar test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

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

	// mat state test
	Mat_Test();

	// opacity test
	Opacity_Test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstOpacity, " << ass.what()
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
	    cout << "**** tstOpacity Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstOpacity on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstOpacity.cc
//---------------------------------------------------------------------------//
