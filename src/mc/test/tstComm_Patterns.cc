//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstComm_Patterns.cc
 * \author Thomas M. Evans
 * \date   Wed May  3 14:36:31 2000
 * \brief  Comm_Patterns test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "DD_Mesh.hh"
#include "../Comm_Patterns.hh"
#include "../General_Topology.hh"
#include "../Rep_Topology.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <string>

using namespace std;

using rtt_mc::Comm_Patterns;
using rtt_mc::Topology;
using rtt_mc::General_Topology;
using rtt_mc::Rep_Topology;
using rtt_dsxx::SP;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__, __FILE__);

//---------------------------------------------------------------------------//
// Comm_Patterns DD Test on 9 cell mesh
//---------------------------------------------------------------------------//

void DD_Comm_Patterns()
{
    if (C4::nodes() != 4)
	return;

    // build a DD topology for a nine cell mesh
    SP<Topology> topology = rtt_mc_test::build_Topology();
    if (topology->get_parallel_scheme() != "DD") ITFAILS;

    // make a cp object
    Comm_Patterns cp;
    
    // nothing should be set yet
    if (cp) ITFAILS;

    // do three permutations of calculating the comm patterns to make sure
    // replication through cycles doesn't F@*% things up
    for (int run = 1; run <= 3; run++)
    {
	// calculate comm_patterns
	cp.calc_patterns(topology);
	if (!cp) ITFAILS;

	// receive comm_patterns

	Comm_Patterns::const_iterator ritor = cp.get_recv_begin();

	// do some checks
	if (C4::node() == 0)
	{
	    if (cp.get_num_recv_procs() != 2) ITFAILS;
	
	    // go through map

	    if (ritor->first != 1)         ITFAILS;
	    if (ritor->second.size() != 2) ITFAILS; 
	    if (ritor->second[0] != 3)     ITFAILS;
	    if (ritor->second[1] != 4)     ITFAILS;
	    ritor++;

	    if (ritor->first != 2)         ITFAILS;
	    if (ritor->second.size() != 1) ITFAILS; 
	    if (ritor->second[0] != 5)     ITFAILS;
	    ritor++;
	}
	else if (C4::node() == 1)
	{
	    if (cp.get_num_recv_procs() != 3) ITFAILS;
	
	    // go through map

	    if (ritor->first != 0)         ITFAILS;
	    if (ritor->second.size() != 2) ITFAILS; 
	    if (ritor->second[0] != 1)     ITFAILS;
	    if (ritor->second[1] != 2)     ITFAILS;
	    ritor++;

	    if (ritor->first != 2)         ITFAILS;
	    if (ritor->second.size() != 2) ITFAILS; 
	    if (ritor->second[0] != 5)     ITFAILS;
	    if (ritor->second[1] != 6)     ITFAILS;
	    ritor++;

	    if (ritor->first != 3)         ITFAILS;
	    if (ritor->second.size() != 1) ITFAILS;  
	    if (ritor->second[0] != 7)     ITFAILS;
	    ritor++;
	}
	else if (C4::node() == 2)
	{
	    if (cp.get_num_recv_procs() != 3) ITFAILS;
	
	    // go through map

	    if (ritor->first != 0)         ITFAILS;
	    if (ritor->second.size() != 1) ITFAILS;
	    if (ritor->second[0] != 2)     ITFAILS;
	    ritor++;

	    if (ritor->first != 1)         ITFAILS;
	    if (ritor->second.size() != 2) ITFAILS; 
	    if (ritor->second[0] != 3)     ITFAILS;
	    if (ritor->second[1] != 4)     ITFAILS;
	    ritor++;

	    if (ritor->first != 3)         ITFAILS;
	    if (ritor->second.size() != 2) ITFAILS;  
	    if (ritor->second[0] != 8)     ITFAILS;  
	    if (ritor->second[1] != 9)     ITFAILS;
	    ritor++;
	}
	else if (C4::node() == 3)
	{
	    if (cp.get_num_recv_procs() != 2) ITFAILS;
	
	    // go through map

	    if (ritor->first != 1)         ITFAILS;
	    if (ritor->second.size() != 1) ITFAILS;
	    if (ritor->second[0] != 4)     ITFAILS;
	    ritor++;

	    if (ritor->first != 2)         ITFAILS;
	    if (ritor->second.size() != 2) ITFAILS; 
	    if (ritor->second[0] != 5)     ITFAILS;
	    if (ritor->second[1] != 6)     ITFAILS;
	    ritor++;
	}

	if (ritor != cp.get_recv_end()) ITFAILS;

	// send comm patterns

	Comm_Patterns::const_iterator sitor = cp.get_send_begin();

	// do some checks
	if (C4::node() == 0)
	{
	    if (cp.get_num_send_procs() != 2) ITFAILS;
	
	    // go through map

	    if (sitor->first != 1)         ITFAILS;
	    if (sitor->second.size() != 2) ITFAILS; 
	    if (sitor->second[0] != 1)     ITFAILS;
	    if (sitor->second[1] != 2)     ITFAILS;
	    sitor++;

	    if (sitor->first != 2)         ITFAILS;
	    if (sitor->second.size() != 1) ITFAILS; 
	    if (sitor->second[0] != 2)     ITFAILS;
	    sitor++;
	}
	else if (C4::node() == 1)
	{
	    if (cp.get_num_send_procs() != 3) ITFAILS;
	
	    // go through map

	    if (sitor->first != 0)         ITFAILS;
	    if (sitor->second.size() != 2) ITFAILS; 
	    if (sitor->second[0] != 3)     ITFAILS;
	    if (sitor->second[1] != 4)     ITFAILS;
	    sitor++;

	    if (sitor->first != 2)         ITFAILS;
	    if (sitor->second.size() != 2) ITFAILS; 
	    if (sitor->second[0] != 3)     ITFAILS;
	    if (sitor->second[1] != 4)     ITFAILS;
	    sitor++;

	    if (sitor->first != 3)         ITFAILS;
	    if (sitor->second.size() != 1) ITFAILS;  
	    if (sitor->second[0] != 4)     ITFAILS;
	    sitor++;
	}
	else if (C4::node() == 2)
	{
	    if (cp.get_num_send_procs() != 3) ITFAILS;
	
	    // go through map

	    if (sitor->first != 0)         ITFAILS;
	    if (sitor->second.size() != 1) ITFAILS;
	    if (sitor->second[0] != 5)     ITFAILS;
	    sitor++;

	    if (sitor->first != 1)         ITFAILS;
	    if (sitor->second.size() != 2) ITFAILS; 
	    if (sitor->second[0] != 5)     ITFAILS;
	    if (sitor->second[1] != 6)     ITFAILS;
	    sitor++;

	    if (sitor->first != 3)         ITFAILS;
	    if (sitor->second.size() != 2) ITFAILS;  
	    if (sitor->second[0] != 5)     ITFAILS;  
	    if (sitor->second[1] != 6)     ITFAILS;
	    sitor++;
	}
	else if (C4::node() == 3)
	{
	    if (cp.get_num_send_procs() != 2) ITFAILS;
	
	    // go through map

	    if (sitor->first != 1)         ITFAILS;
	    if (sitor->second.size() != 1) ITFAILS;
	    if (sitor->second[0] != 7)     ITFAILS;
	    sitor++;

	    if (sitor->first != 2)         ITFAILS;
	    if (sitor->second.size() != 2) ITFAILS; 
	    if (sitor->second[0] != 8)     ITFAILS;
	    if (sitor->second[1] != 9)     ITFAILS;
	    sitor++;
	}

	if (sitor != cp.get_send_end()) ITFAILS;
    }
}

//---------------------------------------------------------------------------//
// Comm_Patterns full replication test.
//---------------------------------------------------------------------------//

void rep_Comm_Patterns()
{
    // make a full replication topology for nine cell mesh
    SP<Topology> topology(new Rep_Topology(9));

    // make a comm_pattern object
    Comm_Patterns cp;

    if (cp) ITFAILS;

    // now build comm patterns.  This should be a null operation because we
    // are in a full replication topology
    cp.calc_patterns(topology);

    if (cp) ITFAILS;
}

//---------------------------------------------------------------------------//
// MAIN
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);
    
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (!C4::node())
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	rep_Comm_Patterns();
	DD_Comm_Patterns();
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	cout << "Dumb ass you screwed up; assertion: " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*****************************************" 
	 << endl;
    if (passed) 
    {
        cout << "**** Comm_Patterns Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "*****************************************" 
	 << endl;
    cout << endl;

    cout << "Done testing Comm_Patterns on node " << C4::node() 
	 << endl;

    C4::Finalize();
}


//---------------------------------------------------------------------------//
//                              end of tstComm_Patterns.cc
//---------------------------------------------------------------------------//
