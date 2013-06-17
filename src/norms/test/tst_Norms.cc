//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/test/tst_Norms.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 09:12:16 2005
 * \brief  Tests Norms.
 * \note   Copyright (C) 2004-2013 Los Alamos National Security, LLC. 
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "norms_test.hh"

#include "../Norms.hh"
#include "../Norms_Labeled.hh"
#include "../Norms_Proc.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_norms::Norms;
using rtt_norms::Norms_Labeled;
using rtt_norms::Norms_Proc;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// Test Norms
void test_Norms()
{
    Norms norms;

    // All norms should be the same if we add just a single quantity, in this
    // case, the processor number + 1.

    size_t x   = rtt_c4::node();
    double xp1 = x + 1.0;

    norms.add(xp1, x);

    UNIT_TEST(soft_equiv(norms.L1(), xp1));
    UNIT_TEST(soft_equiv(norms.L2(), xp1));
    UNIT_TEST(soft_equiv(norms.Linf(), xp1));
    UNIT_TEST(norms.index_Linf() == x);

    // Check accumulating the results back to the host proc.

    norms.comm(0);

    if ( x == 0 )
    {
	size_t num_nodes = rtt_c4::nodes();
	
	UNIT_TEST(norms.index_Linf() == num_nodes-1);
	UNIT_TEST(soft_equiv(norms.Linf(), double(num_nodes)));

	double tL1 = 0.0;
	double tL2 = 0.0;

	for ( size_t i = 0; i < num_nodes; ++i )
	{
	    double xt = i + 1.0;
	    tL1 += xt;
	    tL2 += xt * xt;
	}

	tL1 /= num_nodes;
	tL2 = std::sqrt(tL2 / num_nodes);

	UNIT_TEST(soft_equiv(norms.L1(), tL1));
	UNIT_TEST(soft_equiv(norms.L2(), tL2));
    }

    { // test assignment

	Norms nc;
	nc = norms;
	UNIT_TEST(nc == norms);
    }

    { // test copy ctor.

	Norms nc(norms);
	UNIT_TEST(nc == norms);
    }

    // Done testing

    if ( rtt_norms_test::passed )
    {
        PASSMSG("test_Norms() ok.");
    }
}

// Test Norms_Labeled
void test_Norms_Labeled()
{
    Norms_Labeled norms;

    // All norms should be the same if we add just a single quantity, in this
    // case, the processor number + 1.

    size_t x   = rtt_c4::node();
    double xp1 = x + 1.0;

    Norms_Labeled::Index indx(x);
    ostringstream o;
    o << "proc " << x;
    indx.label = o.str();

    norms.add(xp1, indx);

    UNIT_TEST(soft_equiv(norms.L1(), xp1));
    UNIT_TEST(soft_equiv(norms.L2(), xp1));
    UNIT_TEST(soft_equiv(norms.Linf(), xp1));
    UNIT_TEST(norms.index_Linf() == indx);

    // Check accumulating the results back to the host proc.

    norms.comm(0);

    if ( x == 0 )
    {
	size_t num_nodes = rtt_c4::nodes();
	
	//UNIT_TEST(norms.index_Linf() == num_nodes-1);
	UNIT_TEST(soft_equiv(norms.Linf(), double(num_nodes)));

	double tL1 = 0.0;
	double tL2 = 0.0;

	for ( size_t i = 0; i < num_nodes; ++i )
	{
	    double xt = i + 1.0;
	    tL1 += xt;
	    tL2 += xt * xt;
	}

	tL1 /= num_nodes;
	tL2 = std::sqrt(tL2 / num_nodes);

	UNIT_TEST(soft_equiv(norms.L1(), tL1));
	UNIT_TEST(soft_equiv(norms.L2(), tL2));
    }

    { // test assignment

	Norms_Labeled nc;
	nc = norms;
	UNIT_TEST(nc == norms);
    }

    { // test copy ctor.

	Norms_Labeled nc(norms);
	UNIT_TEST(nc == norms);
    }

    // Done testing

    if ( rtt_norms_test::passed )
    {
        PASSMSG("test_Norms_Labeled() ok.");
    }
}

// Test Norms_Proc
void test_Norms_Proc()
{
    Norms_Proc norms;

    // All norms should be the same if we add just a single quantity, in this
    // case, the processor number + 1.

    size_t x   = rtt_c4::node();
    double xp1 = x + 1.0;

    Norms_Proc::Index indx(x);

    norms.add(xp1, indx);

    UNIT_TEST(soft_equiv(norms.L1(), xp1));
    UNIT_TEST(soft_equiv(norms.L2(), xp1));
    UNIT_TEST(soft_equiv(norms.Linf(), xp1));
    UNIT_TEST(norms.index_Linf() == indx);

    // Check accumulating the results back to the host proc.

    norms.comm(0);

    if ( x == 0 )
    {
	size_t num_nodes = rtt_c4::nodes();
	
	//UNIT_TEST(norms.index_Linf() == num_nodes-1);
	UNIT_TEST(soft_equiv(norms.Linf(), double(num_nodes)));

	double tL1 = 0.0;
	double tL2 = 0.0;

	for ( size_t i = 0; i < num_nodes; ++i )
	{
	    double xt = i + 1.0;
	    tL1 += xt;
	    tL2 += xt * xt;
	}

	tL1 /= num_nodes;
	tL2 = std::sqrt(tL2 / num_nodes);

	UNIT_TEST(soft_equiv(norms.L1(), tL1));
	UNIT_TEST(soft_equiv(norms.L2(), tL2));
    }

    { // test assignment

	Norms_Proc nc;
	nc = norms;
	UNIT_TEST(nc == norms);
    }

    { // test copy ctor.

	Norms_Proc nc(norms);
	UNIT_TEST(nc == norms);
    }

    // Done testing

    if ( rtt_norms_test::passed )
    {
        PASSMSG("test_Norms_Proc() ok.");
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

	test_Norms();
	test_Norms_Labeled();
	test_Norms_Proc();
    }
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing tst_Norms, " 
		  << err.what()
		  << std::endl;
	rtt_c4::finalize();
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tst_Norms, " 
		  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	std::cout << std::endl;
	std::cout <<     "*********************************************" 
                  << std::endl;
	if (rtt_norms_test::passed) 
	{
	    std::cout << "**** tst_Norms Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
	}
    else
	{
	    std::cout << "**** tst_Norms Test: FAILED on " 
                      << rtt_c4::node() 
                      << std::endl;
	}
	std::cout <<     "*********************************************" 
                  << std::endl;
	std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    std::cout << "Done testing tst_Norms on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
// end of tst_Norms.cc
//---------------------------------------------------------------------------//
