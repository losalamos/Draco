//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstRW_Sub_Tally.cc
 * \author Thomas M. Evans
 * \date   Thu Aug 21 16:01:19 2003
 * \brief  Test the Random_Walk_Sub_Tally class.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "../Release.hh"
#include "../Random_Walk_Sub_Tally.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_imc::Random_Walk_Sub_Tally;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void sub_tally_test()
{
    // make sub tally
    Random_Walk_Sub_Tally rwsub;

    if (rwsub.get_accum_n_random_walks() != 0) ITFAILS;
    if (rwsub.get_accum_n_spheres() != 0)      ITFAILS;
    if (rwsub.get_accum_sphere_radii() != 0.0) ITFAILS;
    if (rwsub.get_accum_step_lengths() != 0.0) ITFAILS;

    // do some tallies
    rwsub.accum_n_random_walks();
    rwsub.accum_n_random_walks(4);

    rwsub.accum_sphere_radii(0.1);
    rwsub.accum_sphere_radii(0.4);

    rwsub.accum_step_length(0.8);
    rwsub.accum_step_length(0.7);

    if (rwsub.get_accum_n_random_walks() != 5) ITFAILS;
    if (rwsub.get_accum_n_spheres() != 2)      ITFAILS;
    if (rwsub.get_accum_sphere_radii() != 0.5) ITFAILS;
    if (rwsub.get_accum_step_lengths() != 1.5) ITFAILS;
    
    if (rtt_imc_test::passed)
	PASSMSG("Random_Walk_Sub_Tally ok.");
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
	if (rtt_c4::node() == 0)
	    sub_tally_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstRW_Sub_Tally, " << ass.what()
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
	    cout << "**** tstRW_Sub_Tally Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstRW_Sub_Tally on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstRW_Sub_Tally.cc
//---------------------------------------------------------------------------//
