//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstComm_Dup.cc
 * \author Thomas M. Evans
 * \date   Thu Jul 18 11:10:10 2002
 * \brief  test Communicator Duplication
 * \note   Copyright (C) 2006-2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4_test.hh"
#include "ds++/Release.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_mpi_comm_dup()
{
    // we only run this particular test when mpi is on
#ifdef C4_MPI
    
    Require (rtt_c4::nodes() == 4);

    int node  = rtt_c4::node();
    int snode = 0;

    // split up nodes (two communicators) 0 -> 0, 2 -> 1 and
    // 1 -> 0, 3 -> 1
    MPI_Comm new_comm;
	
    if (node == 1)
    {
	MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &new_comm);
	MPI_Comm_rank(new_comm, &snode);
    }
    else if (node == 3)
    {
	MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &new_comm);
	MPI_Comm_rank(new_comm, &snode);
    }
    else
    {
	MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &new_comm);
	MPI_Comm_rank(new_comm, &snode);
    }
    
    // we haven't set the communicator yet so we should still have 4 nodes
    if (rtt_c4::nodes() != 4) ITFAILS;

    // now dup the communicator on each processor
    rtt_c4::inherit(new_comm);
    
    // each processor should see two nodes
    if (rtt_c4::nodes() != 2) ITFAILS;

    // test data send/receive
    int data = 0;

    // do some tests on each processor
    if (node == 0)
    {
	if (rtt_c4::node() != 0) ITFAILS;

	// set data to 10 and send it out
	data = 10;
	rtt_c4::send(&data, 1, 1, 100);
    }
    else if (node == 1)
    {
	if (rtt_c4::node() != 0) ITFAILS;

	// set data to 20 and send it out
	data = 20;
	rtt_c4::send(&data, 1, 1, 100);
    }
    else if (node == 2)
    {
	if (rtt_c4::node() != 1) ITFAILS;

	if (data != 0) ITFAILS;
	rtt_c4::receive(&data, 1, 0, 100);
	if (data != 10) ITFAILS;
    }
    else if (node == 3) 
    {
	if (rtt_c4::node() != 1) ITFAILS;

	if (data != 0) ITFAILS;
	rtt_c4::receive(&data, 1, 0, 100);
	if (data != 20) ITFAILS;
    }

    // now free the communicator on each processor
    rtt_c4::free_inherited_comm();

    // the free should have set back to COMM_WORLD
    if (rtt_c4::nodes() != 4) ITFAILS;

    {
	rtt_c4::HTSyncSpinLock slock;
	for (int i = 0; i < 10000; i++)
	{
	    continue;
	}

	if (rtt_c4_test::passed)
	{
	    ostringstream m;
	    m << "Communicator duplicated successfully on " << rtt_c4::node();
	    PASSMSG(m.str());
	}
    }

    rtt_c4::global_barrier();

    MPI_Comm_free(&new_comm);

#endif
}

//---------------------------------------------------------------------------//

void test_comm_dup()
{
    // we only run this test scalar
#ifdef C4_SCALAR

    int node = rtt_c4::node();

    // now dup the communicator on each processor
    rtt_c4::inherit(node);

    // check the number of nodes
    if (rtt_c4::nodes() != 1) ITFAILS;

    rtt_c4::free_inherited_comm();

    // check the number of nodes
    if (rtt_c4::nodes() != 1) ITFAILS;

    if (rtt_c4_test::passed)
	PASSMSG("Scalar Comm duplication/free works ok.");

#endif

    // check duping/freeing MPI_COMM_WORLD
#ifdef C4_MPI

    int nodes = rtt_c4::nodes();

    MPI_Comm comm_world = MPI_COMM_WORLD;
    rtt_c4::inherit(comm_world);

    if (rtt_c4::nodes() != nodes) ITFAILS;

    // try a global sum to check
    int x = 10;
    rtt_c4::global_sum(x);
    if (x != 10 * nodes) ITFAILS;

    rtt_c4::free_inherited_comm();

    // we should be back to comm world
    if (rtt_c4::nodes() != nodes) ITFAILS;

    // try a global sum to check
    int y = 20;
    rtt_c4::global_sum(y);
    if (y != 20 * nodes) ITFAILS;

    if (rtt_c4_test::passed)
	PASSMSG("MPI_COMM_WORLD Comm duplication/free works ok."); 
    
#endif
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
		cout << argv[0] << ": version " << rtt_dsxx::release()
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	if (rtt_c4::nodes() == 4)
	    test_mpi_comm_dup();

	test_comm_dup();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstComm_Dup, " << ass.what()
	     << endl;
	rtt_c4::abort();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_c4_test::passed) 
	{
	    cout << "**** tstComm_Dup Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstComm_Dup on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstComm_Dup.cc
//---------------------------------------------------------------------------//
