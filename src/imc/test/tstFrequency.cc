//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstFrequency.cc
 * \author Todd J. Urbatsch
 * \date   Mon Jan 14 11:51:42 2002
 * \brief  Test Frequency class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Frequency.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

using rtt_imc::Gray_Frequency;
using rtt_imc::Multigroup_Frequency;
using rtt_dsxx::soft_equiv;
using rtt_dsxx::SP;

using namespace std;


//---------------------------------------------------------------------------//
// TEST FREQUENCY CLASS
//---------------------------------------------------------------------------//
void tstFreq()
{
    // >>> Test a multigroup Frequency.

    // Make a local vector of group boundaries.
    vector<double> bnds(4, 0.0);
    bnds[0] = 0.1;
    bnds[1] = 0.2;
    bnds[2] = 0.5;
    bnds[3] = 1.0;

    // Make a smart pointer to a Frequency object.
    SP<Multigroup_Frequency> frequency(new Multigroup_Frequency(bnds));

    // This should be a multigroup Frequency.
    if ( frequency->is_gray())          ITFAILS;
    if (!frequency->is_multigroup())    ITFAILS;

    // Test group structure.
    if (frequency->get_num_groups() != bnds.size()-1)         ITFAILS;
    if (frequency->get_num_group_boundaries() != bnds.size()) ITFAILS;

    // Test accessing the entire vector of group boundaries.
    vector<double> accessed_bnds = frequency->get_group_boundaries();
    if (accessed_bnds.size() != frequency->get_num_group_boundaries()) ITFAILS; 

    for (int i = 0; i < frequency->get_num_group_boundaries(); i++)
	if (!soft_equiv(bnds[i], accessed_bnds[i]))                    ITFAILS; 

    // Test getting individual group boundaries.
    for (int g = 1; g <= frequency->get_num_groups(); g++)
    {
	if (!soft_equiv( bnds[g-1], 
			 frequency->get_group_boundaries(g).first )) ITFAILS;
	if (!soft_equiv( bnds[g], 
			 frequency->get_group_boundaries(g).second)) ITFAILS;
    }


    // >>> Test a gray Frequency.

    // Declaration of null local vector of group boundaries.
    vector<double> null_bnds;

    // Make a new gray Frequency.
    SP<Gray_Frequency> gray_frequency(new Gray_Frequency(null_bnds));

    // This should be a gray Frequency.
    if (!gray_frequency->is_gray())          ITFAILS;
    if ( gray_frequency->is_multigroup())    ITFAILS;
}

//---------------------------------------------------------------------------//
// TEST OPERATIONS OF THE FREQUENCY CLASS
//---------------------------------------------------------------------------//
void tstFreq_Ops()
{
    // Useful typedefs.
    typedef rtt_mc::OS_Mesh MT;

    // Make a local vector of group boundaries.
    vector<double> bnds(4, 0.0);
    bnds[0] = 0.1;
    bnds[1] = 0.2;
    bnds[2] = 0.5;
    bnds[3] = 1.0;

    // Make a smart pointer to a new Frequency object.
    SP<Multigroup_Frequency> frequency(new Multigroup_Frequency(bnds));

    // >>> Test the group-search capability.

    // Test the group-search capability.
    if (frequency->find_group_given_a_freq(0.15) != 1)  ITFAILS;
    if (frequency->find_group_given_a_freq(0.25) != 2)  ITFAILS;
    if (frequency->find_group_given_a_freq(0.49) != 2)  ITFAILS;
    if (frequency->find_group_given_a_freq(0.51) != 3)  ITFAILS;

    // Test boundaries.
    if (frequency->find_group_given_a_freq(0.1) != 1)   ITFAILS;
    if (frequency->find_group_given_a_freq(1.0) != 3)   ITFAILS;

    // Test out-of-bounds.
    if (frequency->find_group_given_a_freq(0.05) != -1) ITFAILS;
    if (frequency->find_group_given_a_freq(100.) != -1) ITFAILS;
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

	// test the Frequency class
	tstFreq();

	tstFreq_Ops();

    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstFrequency, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstFrequency Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstFrequency on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                              end of tstFrequency.cc
//---------------------------------------------------------------------------//
