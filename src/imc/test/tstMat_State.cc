//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstMat_State.cc
 * \author Thomas M. Evans
 * \date   Wed Aug 13 10:45:48 2003
 * \brief  Test of Mat_State class.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "imc_test.hh"
#include "../Release.hh"
#include "../Mat_State.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_imc_test::Parser;
using rtt_imc::Mat_State;
using rtt_mc::OS_Builder;
using rtt_mc::OS_Mesh;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::OS_Mesh   MT;
typedef MT::CCSF<double>  ccsf;

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
	PASSMSG("Built simple 2D mesh for testing mat state.");

    return mesh;
}

//---------------------------------------------------------------------------//

SP<Mat_State<MT> > make_mat(SP<MT> mesh)
{
    SP<Mat_State<MT> > mat;
    
    ccsf density(mesh);
    ccsf temp(mesh);
    ccsf cv(mesh);

    density(1) = 3.0;
    density(2) = 3.0;
    density(3) = 3.0;
    density(4) = 4.0;
    density(5) = 5.0;
    density(6) = 6.0;

    temp(1) = 1.0;
    temp(2) = 1.0;
    temp(3) = 1.0;
    temp(4) = 2.0;
    temp(5) = 3.0;
    temp(6) = 5.0;

    cv(1) = 0.1;
    cv(2) = 0.1;
    cv(3) = 0.2;
    cv(4) = 0.1;
    cv(5) = 0.3;
    cv(6) = 0.1;

    mat = new Mat_State<MT>(density, temp, cv);
    return mat;
}

//---------------------------------------------------------------------------//

void mat_state_test()
{
    SP<MT> mesh = build_mesh();

    // make a mat state
    SP<Mat_State<MT> > mat = make_mat(mesh);

    // test it
    
    // density
    {
	if (!soft_equiv(mat->get_rho(1), 3.0)) ITFAILS;
	if (!soft_equiv(mat->get_rho(2), 3.0)) ITFAILS;
	if (!soft_equiv(mat->get_rho(3), 3.0)) ITFAILS;
	if (!soft_equiv(mat->get_rho(4), 4.0)) ITFAILS;
	if (!soft_equiv(mat->get_rho(5), 5.0)) ITFAILS;
	if (!soft_equiv(mat->get_rho(6), 6.0)) ITFAILS;
    }

    // temp
    {
	if (!soft_equiv(mat->get_T(1), 1.0)) ITFAILS
	if (!soft_equiv(mat->get_T(2), 1.0)) ITFAILS;
	if (!soft_equiv(mat->get_T(3), 1.0)) ITFAILS;
	if (!soft_equiv(mat->get_T(4), 2.0)) ITFAILS;
	if (!soft_equiv(mat->get_T(5), 3.0)) ITFAILS;
	if (!soft_equiv(mat->get_T(6), 5.0)) ITFAILS;;
    }

    // spec heat
    {
	if (!soft_equiv(mat->get_spec_heat(1), 0.1)) ITFAILS;
	if (!soft_equiv(mat->get_spec_heat(2), 0.1)) ITFAILS;
	if (!soft_equiv(mat->get_spec_heat(3), 0.2)) ITFAILS;
	if (!soft_equiv(mat->get_spec_heat(4), 0.1)) ITFAILS;
	if (!soft_equiv(mat->get_spec_heat(5), 0.3)) ITFAILS;
	if (!soft_equiv(mat->get_spec_heat(6), 0.1)) ITFAILS;
    }

    if (rtt_imc_test::passed)
    {
	PASSMSG("Mat_State tests ok.");
    }
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
	mat_state_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstMat_State, " << ass.what()
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
	    cout << "**** tstMat_State Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstMat_State on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstMat_State.cc
//---------------------------------------------------------------------------//
