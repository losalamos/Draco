//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstSource_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Dec  8 16:39:33 1999
 * \brief  Source_Builder test file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Source_Builder.hh"
#include "../Release.hh"
#include "mc/OS_Mesh.hh"
#include "c4/global.hh"

#include <vector>
#include <string>
#include <iostream>

using namespace std;

using rtt_imc_test::IMC_Interface;
using rtt_imc::Source_Builder;
using rtt_mc::OS_Mesh;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

//---------------------------------------------------------------------------//
// testing global check equivalences

template<class T>
void test_equivalence(const T value, const T mod_value)
{
    // set a local_value
    T local_value = value;

    // make a source builder
    Source_Builder<OS_Mesh> sb;

    // check if more than 1 node
    if (C4::nodes() > 1)
    {
	// at this point all processors should have the same value
	if (!sb.check_global_equiv(local_value)) ITFAILS;
	

	// now change the first processor's value
	if (C4::node() == 0)
	    local_value = mod_value;

	if (C4::node() > 0)
	{
	    if (!sb.check_global_equiv(local_value)) ITFAILS;
	}
	else
	{
	    if (sb.check_global_equiv(local_value)) ITFAILS;
	}

	// reset all to the same value
	local_value = value;
	if (!sb.check_global_equiv(local_value)) ITFAILS;
	
	// now change the last processor's value
	if (C4::node() == C4::nodes() - 1)
	    local_value = mod_value;

	if (C4::node() == C4::nodes() - 2)
	{
	    if (sb.check_global_equiv(local_value)) ITFAILS;
	}
	else
	{
	    if (!sb.check_global_equiv(local_value)) ITFAILS;
	}
    }
	 
    // reset all to the same value
    local_value = value;
    if (!sb.check_global_equiv(local_value)) ITFAILS;

    // check if more than 2 nodes
    if (C4::nodes() > 2)
    {
	// now change a middle value
	if (C4::node() == C4::nodes()/2)
	    local_value = mod_value;
	
	if (C4::node() == C4::nodes()/2 - 1)
	{
	    if (sb.check_global_equiv(local_value)) ITFAILS;
	}
	else if (C4::node() == C4::nodes()/2)
	{
	    if (sb.check_global_equiv(local_value)) ITFAILS; 
	}
	else
	{
	    if (!sb.check_global_equiv(local_value)) ITFAILS;
	}
    }
	 
    // reset all to the same value
    local_value = value;
    if (!sb.check_global_equiv(local_value)) ITFAILS;

    // check if 1 node (serial)
    if (C4::nodes() == 1)
    {
	// do a serial problem test-->this is trivial but we want to check it
	// anyway
	local_value = mod_value;
	if (!sb.check_global_equiv(local_value)) ITFAILS;
    }
}

void timing()
{
    // set a value
    int value = 10;

    double begin;
    double end;
	
    // make a source builder
    Source_Builder<OS_Mesh> sb;

    if (C4::node() == 0)
	begin = C4::Wtime();

    for (int i = 0; i < 10; i++)
	sb.check_global_equiv(value);

    if (C4::node() == 0)
	end = C4::Wtime();

    if (C4::node() == 0)
	cout << "Ran for " << end-begin << " seconds on " << C4::nodes()
	     << endl;
}

//---------------------------------------------------------------------------//
// main

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl; 
	    C4::Finalize();
	    return 0;
	}

    // test global equivalences
    test_equivalence(10, 11);           // int
    test_equivalence(10.0001, 11.0001); // double
    test_equivalence(10.0001, 10.0002); // double

    // status of test
    cout << endl;
    cout <<     "******************************************" << endl; 
    if (passed) 
    {
        cout << "**** Source_Builder Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "******************************************" << endl;
    cout << endl;

    cout << "Done testing Source_Builder on node: " << C4::node() << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstSource_Builder.cc
//---------------------------------------------------------------------------//
