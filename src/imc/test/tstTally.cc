//----------------------------------*-C++-*----------------------------------//
// tstTally.cc
// Thomas M. Evans
// Wed Apr 28 13:09:29 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Test of Tally class operations
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Release.hh"
#include "../Tally.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "c4/global.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

using rtt_imc_test::Parser;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_imc::Tally;
using rtt_dsxx::SP;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

//---------------------------------------------------------------------------//

void Tally_Test()
{
    // make a mesh and tally
    SP<Parser> parser(new Parser("OS_Input"));
    OS_Builder mb(parser);
    SP<OS_Mesh> mesh = mb.build_Mesh();
    Tally<OS_Mesh> t(mesh);
    
    if (t.num_cells() != 6) ITFAILS;

    // add some stuff and check the tally
    for (int j = 1; j <= 2; j++)
	for (int i = 1; i <= mesh->num_cells(); i++)
	{
	    t.deposit_energy(i, i * 10.0);
	    t.accumulate_ewpl(i, i * 20.0);
	    t.accumulate_cen_info(i, i * 30);
	    t.accum_n_effscat();
	    t.accum_n_thomscat();
	    t.accum_n_killed();
	    t.accum_ew_killed(i * 2.0);
	    t.accum_n_escaped();
	    t.accum_ew_escaped(i * 1.0);
	    t.accum_n_bndcross();
	    t.accum_n_reflections();

	    // checks
	    if (t.get_energy_dep(i) != j * i * 10) ITFAILS;
	    if (t.get_accum_ewpl(i) != j * i * 20) ITFAILS;
	    if (t.get_new_ecen(i) != j * i * 30)   ITFAILS;
	    if (t.get_new_ncen(i) != j * 1)        ITFAILS;
	}
    
    // totals test
    if (t.get_energy_dep_tot() != 20+40+60+80+100+120)  ITFAILS;
    if (t.get_new_ecen_tot() != 60+120+180+240+300+360) ITFAILS;
    if (t.get_new_ncen_tot() != 12)                     ITFAILS;
    if (t.get_ew_escaped() != 2+4+6+8+10+12)            ITFAILS;
    if (t.num_cells() != 6)                             ITFAILS;
}

//---------------------------------------------------------------------------//

void Tally_Test_Evol_Net()
{
    // Make a mesh and field
    SP<Parser> parser(new Parser("OS_Input"));
    OS_Builder mb(parser);
    SP<OS_Mesh> mesh = mb.build_Mesh();
    OS_Mesh::CCSF_double evol_net(mesh);
    for (int i = 1; i <= evol_net.get_Mesh().num_cells(); i++)
	evol_net(i) = i;

    Tally<OS_Mesh> t(mesh, evol_net);
    for (int i = 1; i <= 6; i++)
	if (t.get_evol_net(i) != i) ITFAILS;

    if (t.num_cells() != evol_net.get_Mesh().num_cells()) ITFAILS;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // C4 Init
    C4::Init(argc, argv);

    // this is essentially a scalar test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

    // Tally tests
    Tally_Test();
    Tally_Test_Evol_Net();

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_imc::release() << endl;
	    C4::Finalize();
	    return 0;
	}

    // status of test
    cout << endl;
    cout <<     "*********************************" << endl;
    if (passed) 
    {
        cout << "**** Tally Self Test: PASSED ****" << endl;
    }
    cout <<     "*********************************" << endl;
    cout << endl;

    cout << "Done testing Tally." << endl;

    // C4 end
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                              end of tstTally.cc
//---------------------------------------------------------------------------//
