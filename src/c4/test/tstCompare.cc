//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstCompare.cc
 * \author Mike Buksas
 * \date   Thu May  1 14:47:00 2008
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "c4_test.hh"
#include "ds++/Release.hh"
#include "../global.hh"
#include "../SpinLock.hh"

#include "../Compare.hh"

using namespace std;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

template<class T>
void test_equivalence(const T value, const T alt_value)
{
    T local_value = value;

    // Test requires more than one node:
    if (rtt_c4::nodes() > 1)
    {
	// at this point all processors should have the same value
	if (!check_global_equiv(local_value)) ITFAILS;
	

	// now change the first processor's value
	if (rtt_c4::node() == 0)
	    local_value = alt_value;

	if (rtt_c4::node() > 0)
	{
	    if (!check_global_equiv(local_value)) ITFAILS;
	}
	else
	{
	    if (check_global_equiv(local_value)) ITFAILS;
	}

	// Reset to given value
	local_value = value;
	if (!check_global_equiv(local_value)) ITFAILS;
	
	// Change the last processor:
	if (rtt_c4::node() == rtt_c4::nodes() - 1)
	    local_value = alt_value;

	if (rtt_c4::node() == rtt_c4::nodes() - 2)
	{
	    if (check_global_equiv(local_value)) ITFAILS;
	}
	else
	{
	    if (!check_global_equiv(local_value)) ITFAILS;
	}
    }
	 
    // Reset to given value
    local_value = value;
    if (!check_global_equiv(local_value)) ITFAILS;

    // Test valid on two nodes or more:
    if (rtt_c4::nodes() > 2)
    {
	// Change a middle value
	if (rtt_c4::node() == rtt_c4::nodes()/2)
	    local_value = alt_value;
	
	if (rtt_c4::node() == rtt_c4::nodes()/2 - 1)
	{
	    if (check_global_equiv(local_value)) ITFAILS;
	}
	else if (rtt_c4::node() == rtt_c4::nodes()/2)
	{
	    if (check_global_equiv(local_value)) ITFAILS; 
	}
	else
	{
	    if (!check_global_equiv(local_value)) ITFAILS;
	}
    }
	 
    // Reset
    local_value = value;
    if (!check_global_equiv(local_value)) ITFAILS;

    // Check 1 node. trivial, but check anyway.
    if (rtt_c4::nodes() == 1)
    {
	local_value = alt_value;
	if (!check_global_equiv(local_value)) ITFAILS;
    }

}
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " 
		     << rtt_dsxx::release()
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// test global equivalences
	test_equivalence(10, 11);           // int
	test_equivalence(10.0001, 11.0001); // double
	test_equivalence(10.0001, 10.0002); // double
    test_equivalence(uint64_t(10000000000), uint64_t(200000000000)); // uint64_t

    }
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing tstCompare, " 
		  << err.what()
		  << std::endl;
	rtt_c4::abort();
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tstCompare, " 
		  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
	rtt_c4::abort();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	std::cout << std::endl;
	std::cout <<     "*********************************************" 
                  << std::endl;
	if (rtt_c4_test::passed) 
	{
	    std::cout << "**** tstCompare Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
	}
	std::cout <<     "*********************************************" 
                  << std::endl;
	std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    std::cout << "Done testing tstCompare on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   


//---------------------------------------------------------------------------//
//                        end of tstCompare.cc
//---------------------------------------------------------------------------//
